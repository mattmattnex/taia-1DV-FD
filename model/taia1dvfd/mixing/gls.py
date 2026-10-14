import math
import torch
from taia1dvfd.core import Model, pcr_solve

def gls_init(model: Model):
    cfg   = model.config
    state = model.state

    cfg.gls = type('GLSParams', (), {})
    p = cfg.gls

    _stabf_parameters(model)
    _model_parameters(model)

    # --- Initial fields ---
    l0    = cfg.gls_lsc0
    diss0 = p.cm0 ** 3 * p.tke_min ** 1.5 / l0
    psi0  = p.cm0 ** p.rp * p.tke_min ** p.rm * l0 ** p.rn

    state.az_m[:] = cfg.gls_azm_min
    state.az_t[:] = cfg.gls_azt_min
    state.tke[:]  = p.tke_min
    state.psi[:]  = psi0
    state.lsc[:]  = l0
    state.diss[:] = diss0

def gls_step(u_star, model: Model):
    state = model.state

    # --- Before-step values ---
    tke_old = state.tke.clone()
    psi_old = state.psi.clone()

    # --- Production/dissipation ---
    P_w, M2_E = _compute_production(model)
    diss_old = _compute_dissipation(tke_old, psi_old, model)

    # --- Solve k, then psi ---
    tke_new = _tke_equation(u_star, tke_old, diss_old, P_w, model)
    psi_new = _psi_equation(u_star, tke_old, psi_old, diss_old, P_w, tke_new, model)

    # --- Diagnostics from the new state ---
    diss_new = _compute_dissipation(tke_new, psi_new, model)
    lsc_new  = _compute_length_scale(tke_new, diss_new, model)

    Akv_new, Akt_new = _compute_mixing(tke_new, diss_new, M2_E, model)

    state.tke[:]  = tke_new
    state.psi[:]  = psi_new
    state.diss[:] = diss_new
    state.lsc[:]  = lsc_new
    state.az_m[:] = Akv_new
    state.az_t[:] = Akt_new


# =====================================================================
#  Internal functions
# =====================================================================
def _compute_production(model):
    cfg, grd, st = model.config, model.grid, model.state
    u_tot  = (st.ue + st.us) if cfg.gls_wave_mix else st.ue
    du     = st.ue[:, :-1] - st.ue[:, 1:]
    dut    = u_tot[:, :-1] - u_tot[:, 1:]
    a_half = 0.5 * (st.az_m[:, :-1] + st.az_m[:, 1:])
    P_w  = grd.to_w(a_half * (du * dut.conj()).real / grd.dz**2).clamp(min=0.0)   # k & psi production
    M2_E = grd.to_w((du * du.conj()).real / grd.dz**2).clamp(min=cfg.gls_M2_min)  # alpha_m only
    return P_w, M2_E

def _compute_dissipation(tke, psi, model: Model):
    """eps = cm0^(3+p/n) k^(1.5+m/n) psi^(-1/n)."""
    p = model.config.gls
    diss = (p.cm0 ** p.e1) * (tke ** p.e2) * (psi ** p.e3)
    return torch.clamp(diss, min=p.diss_min)

def _compute_length_scale(tke, diss, model: Model):
    """l = cm0^3 k^{3/2} / eps."""
    cfg = model.config
    p = cfg.gls
    lsc = (p.cm0 ** 3) * (tke ** 1.5) / diss
    return torch.clamp(lsc, min=cfg.gls_lsc_min)

def _build_operator(sink, sigma, model: Model):
    """I + dt*diag(sink) - (dt/sigma) * d/dz( avm d/dz ), with the face
    diffusivity 0.5*(avm[j]+avm[j+1]) built here (NOT ops.d_azmdz, whose
    boundary rows are overwritten in place by the momentum step and whose
    coefficients lag one time step)."""
    cfg = model.config
    grd = model.grid
    azm = model.state.az_m

    dt, nz  = cfg.dt, cfg.nz
    inv_dz2 = grd.inv_dz2

    I = torch.eye(nz + 1, **cfg.farr_kwargs)
    L = I + torch.diag_embed((dt * sink).squeeze(-1))

    azm_v = azm[:, :, 0]
    jk = torch.arange(1, nz, device=azm.device)

    a_up = 0.5 * (azm_v[:, jk + 1] + azm_v[:, jk    ]) * inv_dz2 / sigma
    a_lw = 0.5 * (azm_v[:, jk    ] + azm_v[:, jk - 1]) * inv_dz2 / sigma

    L[:, jk, jk + 1] += -dt * a_up
    L[:, jk, jk - 1] += -dt * a_lw
    L[:, jk, jk    ] +=  dt * (a_up + a_lw)
    return L

def _set_dirichlet(L, RHS, row, value):
    L[:, row, :]   = 0.0
    L[:, row, row] = 1.0
    RHS[:, row, :] = value

def _set_neumann(L, RHS, row, nbr, flux, dt, dz):
    """Decouple `row` from its boundary neighbour `nbr` and inject `flux`."""
    L[:, row, row] += L[:, row, nbr]
    L[:, row, nbr]  = 0.0
    RHS[:, row, :] += dt * flux / dz

def _tke_equation(u_star, tke_b, diss_b, prod, model: Model):
    cfg = model.config
    p   = cfg.gls
    dt  = cfg.dt
    dz  = model.grid.dz

    N2   = model.state.N2
    az_t = model.state.az_t

    # --- Sources / sinks ---
    P   =  prod                            # shear production
    P_b = -az_t * N2                       # buoyancy term

    B_expl = torch.clamp(P_b, min=0.0)     # source (unstable)  -> explicit RHS
    B_impl = torch.clamp(P_b, max=0.0)     # sink   (stable)    -> implicit

    sink = (diss_b - B_impl) / tke_b       # Patankar: diagonal stays > 1

    # --- System ---
    L   = _build_operator(sink, p.sig_tke, model)
    RHS = tke_b + dt * (P + B_expl)

    # --- Boundary conditions [Warner Eq. 52] ---
    tke_sfc = torch.clamp(p.cm0_inv2 * u_star ** 2, min=p.tke_min)
    tke_bot = torch.full_like(tke_sfc, p.tke_min)      # zero bottom friction

    if cfg.gls_bc == 'dirichlet':          # two rows per boundary (GOTM/NEMO)
        _set_dirichlet(L, RHS,  0, tke_sfc)
        _set_dirichlet(L, RHS,  1, tke_sfc)
        _set_dirichlet(L, RHS, -1, tke_bot)
        _set_dirichlet(L, RHS, -2, tke_bot)
    else:                                  # Neumann: zero k-flux, no wave brk.
        _set_dirichlet(L, RHS,  0, tke_sfc)
        _set_dirichlet(L, RHS, -1, tke_bot)
        _set_neumann(L, RHS,  1,  0, 0.0, dt, dz)
        _set_neumann(L, RHS, -2, -1, 0.0, dt, dz)

    # tke_new = torch.linalg.solve(L, RHS)
    tke_new = pcr_solve(L, RHS)
    return torch.clamp(tke_new, min=p.tke_min)

def _psi_equation(u_star, tke_b, psi_b, diss_b, prod, tke_new, model: Model):
    cfg = model.config
    p   = cfg.gls
    dt  = cfg.dt
    dz  = model.grid.dz

    N2   = model.state.N2
    az_t = model.state.az_t

    # --- Coefficients (before values, NEMO/ROMS convention) ---
    inv_tke = 1.0 / tke_b
    ratio   = psi_b * inv_tke

    # Shear production
    P = p.beta1 * ratio * prod

    # Buoyancy: c3- on the stable part, c3+ on the unstable part
    B = -ratio * az_t * (p.beta3m * torch.clamp(N2, min=0.0)
                       + p.beta3p * torch.clamp(N2, max=0.0))

    B_expl = torch.clamp(B, min=0.0)
    B_impl = torch.clamp(B, max=0.0)

    # --- Wall function (k-kl only) [Warner Eq. 20, both boundaries] ---
    Fwall = _get_Fwall(tke_b, diss_b, model) if cfg.gls_model == 'k-kl' else 1.0

    diss_coeff = p.beta2 * inv_tke * diss_b * Fwall
    buoy_coeff = -B_impl / (psi_b + cfg._eps_sf)
    sink = diss_coeff + buoy_coeff

    # --- System ---
    L   = _build_operator(sink, p.sig_psi, model)
    RHS = psi_b + dt * (P + B_expl)

    # --- Boundary conditions [Warner Eqs. 53-54] ---
    psi_s0, psi_s1, psi_b0, psi_b1, flx_s, flx_b = \
        _get_psi_bc(u_star, tke_new, model)

    if cfg.gls_bc == 'dirichlet':
        _set_dirichlet(L, RHS,  0, psi_s0)
        _set_dirichlet(L, RHS,  1, psi_s1)
        _set_dirichlet(L, RHS, -1, psi_b0)
        _set_dirichlet(L, RHS, -2, psi_b1)
    else:
        _set_dirichlet(L, RHS,  0, psi_s0)
        _set_dirichlet(L, RHS, -1, psi_b0)
        _set_neumann(L, RHS,  1,  0, flx_s, dt, dz)
        _set_neumann(L, RHS, -2, -1, flx_b, dt, dz)

    # psi_new = torch.linalg.solve(L, RHS)
    psi_new = pcr_solve(L, RHS)
    psi_new = torch.clamp(psi_new, min=cfg.gls_psi_min)

    # --- Length-scale limitation ---
    # Geometric bound (law of the wall, from the true surface) everywhere,
    # Galperin (1988) bound under stable stratification.
    depth_w = cfg.Hmin + torch.arange(cfg.nz + 1, **cfg.farr_kwargs).view(1, -1, 1) * dz
    lsc_geo = cfg.vonKar * (depth_w + cfg.z0_surf)        # (1, nz+1, 1), broadcasts over ne
    
    stable  = N2 > 0
    lsc_gal = cfg.gls_galp * torch.sqrt(2.0 * tke_new / torch.clamp(N2, min=cfg.gls_N2_min))
    lsc_gal = torch.where(stable, lsc_gal, torch.full_like(lsc_gal, float('inf')))

    lsc_max = torch.minimum(lsc_geo, lsc_gal)
    psi_lim = (p.cm0 ** p.rp) * (tke_new ** p.rm) * (lsc_max ** p.rn)

    interior = torch.zeros_like(stable); interior[:, 2:-2] = True
    if p.rn < 0:   # psi decreases with l -> floor on psi
        psi_new = torch.where(interior, torch.maximum(psi_new, psi_lim), psi_new)
    else:
        psi_new = torch.where(interior, torch.minimum(psi_new, psi_lim), psi_new)

    return psi_new


def _compute_mixing(tke, diss, M2, model: Model):
    cfg = model.config

    alpha_m, alpha_n = _get_alpha_mn(tke, diss, M2, model)
    c_mu, c_mu_prim  = _get_c_mu(model, alpha_m, alpha_n)

    cff_visc = tke ** 2 / diss
    Akv = torch.clamp(cff_visc * c_mu,      min=cfg.gls_azm_min, max=cfg.az_max)
    Akt = torch.clamp(cff_visc * c_mu_prim, min=cfg.gls_azt_min, max=cfg.az_max)
    return Akv, Akt


# =====================================================================
#  Boundary conditions
# =====================================================================

def _get_psi_bc(u_star, tke, model: Model):
    """Dirichlet values [Warner Eq. 53] at the two rows nearest each boundary
    (distances d = 0 and d = dz on the w-grid) and analytical fluxes
    [Warner Eq. 54] evaluated at d = dz/2 for the Neumann option.

      psi(d)  = cm0^p k^m ( kappa (z0 + d) )^n
      F_psi   = -n cm0^{p+1} kappa / sig_psi  k^{m+1/2} ( kappa (z0 + d) )^n
    """
    cfg = model.config
    p   = cfg.gls
    dz  = model.grid.dz
    d0 = cfg.Hmin

    cm0    = p.cm0
    vonKar = cfg.vonKar
    rp, rm, rn = p.rp, p.rm, p.rn
    cff_flx = -rn * cm0 ** (rp + 1.0) * vonKar / p.sig_psi

    # --- Surface (Charnock roughness) ---
    z0_s = torch.full_like(u_star, cfg.z0_surf)
    k_s0, k_s1 = tke[:, 0, :], tke[:, 1, :]

    psi_s0 = cm0**rp * k_s0**rm * (vonKar * (z0_s + d0     ))**rn
    psi_s1 = cm0**rp * k_s1**rm * (vonKar * (z0_s + d0 + dz))**rn
    lgth_s = vonKar * (z0_s + d0 + 0.5 * dz)
    
    k_sh   = 0.5 * (k_s0 + k_s1)
    flx_s  = cff_flx * (k_sh ** (rm + 0.5)) * (lgth_s ** rn)

    # --- Bottom (no friction: z0 floor only) ---
    z0_b = cfg.gls_zob_min
    k_b0, k_b1 = tke[:, -1, :], tke[:, -2, :]

    psi_b0 = (cm0 ** rp) * (k_b0 ** rm) * (vonKar *  z0_b      ) ** rn
    psi_b1 = (cm0 ** rp) * (k_b1 ** rm) * (vonKar * (z0_b + dz)) ** rn

    lgth_b = vonKar * (z0_b + 0.5 * dz)
    k_bh   = 0.5 * (k_b0 + k_b1)
    flx_b  = cff_flx * (k_bh ** (rm + 0.5)) * (lgth_b ** rn)

    psi_min = cfg.gls_psi_min
    return (torch.clamp(psi_s0, min=psi_min), torch.clamp(psi_s1, min=psi_min),
            torch.clamp(psi_b0, min=psi_min), torch.clamp(psi_b1, min=psi_min),
            flx_s, flx_b)

def _get_Fwall(tke, diss, model: Model):
    """MY wall proximity function, k-kl only.  [Warner Eq. 20]

      F = 1 + E2 (L / kappa)^2 (1/d_s)^2

    with d_s, d_b the distances to the surface and bottom w-points.
    """
    cfg = model.config
    p   = cfg.gls
    grd = model.grid
    nz  = cfg.nz

    dz  = grd.dz
    d_s = cfg.Hmin + torch.arange(nz + 1, **cfg.farr_kwargs).view(1, nz + 1, 1) * dz
    d_s += cfg.z0_surf

    L = torch.clamp(p.cm0 ** 3 * tke ** 1.5 / diss, min= cfg.gls_lsc_min)
    F = 1.0 + 1.33 * (L / cfg.vonKar) ** 2 * (1.0 / d_s) ** 2

    return torch.clamp(F, max=1e4)

# =====================================================================
#  Stability functions
# =====================================================================

def _get_c_mu(model: Model, a_m, a_n):
    """c_mu(a_m, a_n) and c_mu'(a_m, a_n), GOTM weak-equilibrium form."""
    cfg = model.config
    p   = cfg.gls

    den = (p._d0 + p._d1 * a_n + p._d2 * a_m + p._d3 * a_n * a_m
           + p._d4 * a_n ** 2 + p._d5 * a_m ** 2)

    c_mu      = (p._n0  + p._n1  * a_n + p._n2  * a_m) / (den + cfg._eps_sf)
    c_mu_prim = (p._nb0 + p._nb1 * a_n + p._nb2 * a_m) / (den + cfg._eps_sf)

    c_mu      = torch.clamp(c_mu,      min=cfg.gls_cmu_min)
    c_mu_prim = torch.clamp(c_mu_prim, min=cfg.gls_cmu_prim_min)
    return c_mu, c_mu_prim

def _alpha_m_max(model: Model, a_n):
    cfg = model.config
    p   = cfg.gls
    num = p._lam0 + p._lam1 * a_n + p._lam2 * a_n ** 2 + p._lam3 * a_n ** 3
    den = p._lam4 + p._lam5 * a_n + p._lam6 * a_n ** 2
    return torch.clamp(num / (den + cfg._eps_sf), min=0.0, max=cfg._lim_sf)

def _get_alpha_mn(tke, diss, M2, model: Model):
    """alpha_m = (k/eps)^2 M^2,  alpha_n = (k/eps)^2 N^2, with limiters."""
    cfg = model.config
    p   = cfg.gls
    N2  = model.state.N2

    cff = (tke / diss) ** 2.0
    a_m = cff * M2
    a_n = cff * N2

    a_n = torch.clamp(a_n, max=cfg._lim_sf, min=p.alpha_n_min)
    a_m = torch.clamp(a_m, max=_alpha_m_max(model, a_n))
    return a_m, a_n

# =====================================================================
#  Initializing
# =====================================================================

def _stabf_parameters(model: Model):
    cfg = model.config
    c1,  c2,  c3,  c4,  c5,  c6  = cfg.gls_stabf_params[cfg.gls_stabf]["c"]
    cb1, cb2, cb3, cb4, cb5, cb6 = cfg.gls_stabf_params[cfg.gls_stabf]["cb"]

    a1 = 2.0 / 3.0 - 0.5 * c2
    a2 = 1.0       - 0.5 * c3
    a3 = 1.0       - 0.5 * c4
    a5 = 0.5       - 0.5 * c6
    ab1 = 1.0 - cb2
    ab2 = 1.0 - cb3
    ab3 = 2.0 * (1.0 - cb4)
    ab5 = 2.0 * (1.0 - cb5) * cb6

    nn = 0.5 * c1
    nb = cb1

    c_mu_0 = ((a2 ** 2 - 3.0 * a3 ** 2 + 3.0 * a1 * nn) / (3.0 * nn ** 2)) ** 0.25
    cfg.gls.cm0 = c_mu_0
    cfg.gls.cm0_inv2 = 1.0 / (c_mu_0 ** 2)

    # --- Denominator coefficients ---
    cfg.gls._d0 = 36.0 * nn ** 3 * nb ** 2
    cfg.gls._d1 = 84.0 * a5 * ab3 * nn ** 2 * nb + 36.0 * ab5 * nn ** 3 * nb
    cfg.gls._d2 = (9.0 * (ab2 ** 2 - ab1 ** 2) * nn ** 3
                   - 12.0 * (a2 ** 2 - 3.0 * a3 ** 2) * nn * nb ** 2)
    cfg.gls._d3 = (12.0 * a5 * ab3 * (a2 * ab1 - 3.0 * a3 * ab2) * nn
                   + 12.0 * a5 * ab3 * (a3 ** 2 - a2 ** 2) * nb
                   + 12.0 * ab5 * (3.0 * a3 ** 2 - a2 ** 2) * nn * nb)
    cfg.gls._d4 = 48.0 * a5 ** 2 * ab3 ** 2 * nn + 36.0 * a5 * ab3 * ab5 * nn ** 2
    cfg.gls._d5 = 3.0 * (a2 ** 2 - 3.0 * a3 ** 2) * (ab1 ** 2 - ab2 ** 2) * nn

    # --- c_mu numerator coefficients ---
    cfg.gls._n0 = 36.0 * a1 * nn ** 2 * nb ** 2
    cfg.gls._n1 = (-12.0 * a5 * ab3 * (ab1 + ab2) * nn ** 2
                   + 8.0 * a5 * ab3 * (6.0 * a1 - a2 - 3.0 * a3) * nn * nb
                   + 36.0 * a1 * ab5 * nn ** 2 * nb)
    cfg.gls._n2 = 9.0 * a1 * (ab2 ** 2 - ab1 ** 2) * nn ** 2

    # --- c_mu' numerator coefficients ---
    cfg.gls._nb0 = 12.0 * ab3 * nn ** 3 * nb
    cfg.gls._nb1 = 12.0 * a5 * ab3 ** 2 * nn ** 2
    cfg.gls._nb2 = (9.0 * a1 * ab3 * (ab1 - ab2) * nn ** 2
                    + (6.0 * a1 * (a2 - 3.0 * a3)
                       - 4.0 * (a2 ** 2 - 3.0 * a3 ** 2)) * ab3 * nn * nb)

    # --- alpha_m upper-bound coefficients ---
    p = cfg.gls
    cfg.gls._lam0 = p._d0 * p._n0
    cfg.gls._lam1 = p._d0 * p._n1 + p._d1 * p._n0
    cfg.gls._lam2 = p._d1 * p._n1 + p._d4 * p._n0
    cfg.gls._lam3 = p._d4 * p._n1
    cfg.gls._lam4 = p._d2 * p._n0
    cfg.gls._lam5 = p._d2 * p._n1 + p._d3 * p._n0
    cfg.gls._lam6 = p._d3 * p._n1

    # --- alpha_n lower bound ---
    A = p._d4 + p._nb1
    B = p._d1 + p._nb0
    C = p._d0
    disc = B * B - 4.0 * A * C
    root = (-B + math.sqrt(max(disc, 0.0))) / (2.0 * A)
    cfg.gls.alpha_n_min = 0.73 * root


def _model_parameters(model: Model):
    cfg = model.config
    prm = cfg.gls_model_params[cfg.gls_model]

    # Exponents
    cfg.gls.rp = prm['rp']
    cfg.gls.rm = prm['rm']
    cfg.gls.rn = prm['rn']
    cfg.gls.e1 =  3.0 + prm['rp'] / prm['rn']
    cfg.gls.e2 =  1.5 + prm['rm'] / prm['rn']
    cfg.gls.e3 = -1.0 / prm['rn']

    # Coefficients
    cfg.gls.sig_tke = prm['sig_tke']
    cfg.gls.sig_psi = prm['sig_psi']
    cfg.gls.beta1   = prm['beta1']
    cfg.gls.beta2   = prm['beta2']
    cfg.gls.beta3p  = prm['beta3p']
    cfg.gls.beta3m  = prm['beta3m'][cfg.gls_stabf]

    cfg.gls.tke_min  = cfg.gls_tke_min
    cfg.gls.diss_min = cfg.gls_diss_min