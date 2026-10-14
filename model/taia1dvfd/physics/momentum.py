import torch
from taia1dvfd.core import Model, pcr, tridiag_apply, ddz

def momentum_equation_init(model: Model) -> None:
    pass

def momentum_equation_step(model: Model) -> None:
    """
    Linear-implicit Euler-Maruyama step of the LU Ekman-Stokes momentum equation,
    in banded form: only the three diagonals of the operators are built (ne, n),
    never the dense (ne, n, n) matrices.

      d u = [ -if (u + us) + dz( a dz (u + us) ) ] dt  -  [ if sig_h + sig_z dz(u + us) ] dB
    (wave-mixing terms dropped when cfg.wave_mixing is False).
    """
    state = model.state
    cfg   = model.config
    grd   = model.grid

    # --- Variables (ne, n) with n = nz+1 face nodes ---
    ue    = state.ue[..., 0]
    us    = state.us[..., 0]
    ua    = state.ua
    azm   = state.az_m[..., 0]
    sig_z = state.sig_z_m[..., 0]
    sig_h = state.sig_h[..., 0]
    dB    = state.dB_w[..., 0]

    # --- Parameters ---
    dz    = grd.dz
    dt    = cfg.dt
    f     = cfg.f
    cdrho = cfg.cdrho

    # ==================================================================
    # --- Diffusion operator D = dz( a dz . ), Neumann rows (dz/2 volume) ---
    # ==================================================================
    a_half = 0.5 * (azm[:, :-1] + azm[:, 1:])         # (ne, nz): a at i+1/2

    D_up = a_half / dz**2                             # D[i, i+1], rows 0..nz-1
    D_lo = a_half / dz**2                             # D[i+1, i], rows 1..nz
    D_up[:, 0]  = 2.0 * D_up[:, 0]                    # surface node: half control volume
    D_lo[:, -1] = 2.0 * D_lo[:, -1]                   # bottom node:  half control volume

    D_dg = torch.zeros_like(azm)
    D_dg[:, :-1] -= D_up
    D_dg[:, 1:]  -= D_lo

    # --- A = D - i f I   (Coriolis on the diagonal) ---
    A_dg = D_dg - 1j * f

    # ==================================================================
    # --- RHS ---
    # ==================================================================
    if cfg.wave_mixing:
        Fk    = tridiag_apply(D_lo, A_dg, D_up, us)       # A us : Coriolis-Stokes + wave mixing
        u_tot = ue + us
    else:
        Fk    = -1j * f * us                             # Coriolis-Stokes only
        u_tot = ue

    rhs = ue + dt * Fk

    if cfg.vel_stochastic:
        Gk  = 1j * f * sig_h + sig_z * ddz(u_tot, dz)
        rhs = rhs - Gk * dB

    # --- Boundary forcing (fluxes over the dz/2 boundary volumes) ---
    tau_w = (cdrho * torch.abs(ua) * ua)[:, 0, 0]        # (ne,)
    tau_s = azm[:, 0] * us[:, 0] / cfg.hs                # (ne,)

    Fsurf = tau_w + tau_s if cfg.wave_mixing else tau_w
    Fbotm = torch.zeros_like(Fsurf)

    rhs[:,  0] += dt * Fsurf / (0.5 * dz)
    rhs[:, -1] += dt * Fbotm / (0.5 * dz)

    # ==================================================================
    # --- Solve (I - dt A) u^{n+1} = rhs ---
    # ==================================================================
    L_lo = -dt * D_lo
    L_up = -dt * D_up
    L_dg = 1.0 - dt * A_dg

    state.ue[:] = pcr(L_lo, L_dg, L_up, rhs.unsqueeze(-1))