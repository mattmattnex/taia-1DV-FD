import torch
from taia1dvfd.core import Model, Grid, Config

def tke_init(model: Model) -> None:
    state = model.state
    cfg   = model.config

    state.az_m[:]  = cfg.tke_azm_min
    state.az_t[:]  = cfg.tke_azt_min
    state.tke[:]   = cfg.tke_min
    state.dissl[:] = cfg.tke_min ** 0.5 / cfg.tke_rmxl_min

def tke_step(u_star: torch.Tensor, model: Model) -> None:
    state = model.state
    cfg   = model.config
    grd   = model.grid

    en    = state.tke
    azm   = state.az_m
    azt   = state.az_t
    N2    = state.N2
    dissl = state.dissl

    fkw     = cfg.farr_kwargs
    dt      = cfg.dt
    nz      = cfg.nz
    ne      = cfg.ne
    inv_dz2 = grd.inv_dz2

    Ck    = cfg.tke_ediff
    Ceps  = cfg.tke_ediss
    emin  = cfg.tke_min
    emin0 = cfg.tke_surf_min

    # --- Shear ---
    u_tot  = (state.ue + state.us) if cfg.tke_wave_mix else state.ue
    du     = state.ue[:, :-1] - state.ue[:, 1:]               # (ne, nz) differences between face nodes
    dut    = u_tot[:, :-1]    - u_tot[:, 1:]
    a_half = 0.5 * (azm[:, :-1] + azm[:, 1:])                 # same a_{i+1/2} as the momentum operator
    P_c    = a_half * (du * dut.conj()).real / grd.dz**2      # exactly what d_azmdz removes from u
    sh2    = grd.to_w(P_c)                     # production at faces
    
    P_c_E = a_half * (du * du.conj()).real / grd.dz**2
    sh2_E = grd.to_w(P_c_E) # for Richardson number

    # ====================
    # --- TKE equation ---
    # ====================
    Idc = torch.eye(nz + 1, **fkw)
    L   = Idc.unsqueeze(0).repeat(ne, 1, 1)

    azm_v = azm[:, :, 0]
    jk = torch.arange(1, nz, device=azm.device)

    # NEMO clamps the SUM of adjacent avm at 2e-5 (effective floor 1e-5 on the mean)
    a_up = 0.5 * torch.clamp(azm_v[:, jk + 1] + azm_v[:, jk],     min=cfg.tke_azm_floor) * inv_dz2
    a_lw = 0.5 * torch.clamp(azm_v[:, jk]     + azm_v[:, jk - 1], min=cfg.tke_azm_floor) * inv_dz2

    L[:, jk, jk + 1] += -dt * a_up
    L[:, jk, jk - 1] += -dt * a_lw
    
    # Linearization of Ceps e^{3/2}/l :
    L[:, jk, jk]     +=  dt * (a_up + a_lw) + 1.5 * dt * Ceps * dissl[:, jk, 0]

    rhs = en.clone()
    rhs[:, 1:nz, 0] += dt * (sh2[:, 1:nz, 0]
                             - azt[:, 1:nz, 0] * N2[:, 1:nz, 0]
                             + 0.5 * Ceps * dissl[:, 1:nz, 0] * en[:, 1:nz, 0])

    # --- BCs (Dirichlet) ---
    en_surf = torch.clamp(cfg.tke_ebb * u_star ** 2, min=emin0)

    L[:, 0, :]   = 0.0
    L[:, 0, 0]   = 1.0
    rhs[:, 0, :] = en_surf

    L[:, -1, :]   = 0.0
    L[:, -1, -1]  = 1.0
    rhs[:, -1, :] = emin

    en_new = torch.linalg.solve(L, rhs)
    en_new = torch.clamp(en_new, min=emin)
    state.tke[:] = en_new

    # ==========================================================
    # --- mixing lengths from the NEW en  ->  avm, avt, dissl ---
    # ==========================================================
    lk, le = _mixing_lengths(en_new, N2, cfg, grd)
    sq_en   = torch.sqrt(en_new)
    state.dissl[:] = sq_en / le
    
    azm_raw = Ck * lk * sq_en
    azm_new = torch.clamp(azm_raw, min=cfg.tke_azm_min)
    
    # --- Prandtl number ---
    Ri   = N2 * azm_new / (sh2_E + cfg.tke_bshear)
    Ri   = torch.where(N2 <= 0.0, torch.zeros_like(Ri), Ri)
    pdlr = (cfg.tke_ri_cri / Ri.clamp(min=cfg.tke_ri_cri)).clamp(min=0.1)        

    azt_new = (azm_raw * pdlr).clamp(min=cfg.tke_azt_min)
    
    state.az_m[:] = azm_new
    state.az_t[:] = azt_new
        
def _mixing_lengths(en: torch.Tensor, N2: torch.Tensor, cfg: Config, grd: Grid):
    nz       = cfg.nz
    dz       = grd.dz
    rmxl_min = cfg.tke_rmxl_min

    N2c = N2.clamp(min=cfg.tke_N2_min)
    mxl = torch.sqrt(2.0 * en.clamp(min=cfg.tke_min) / N2c)
    mxl = mxl.clamp(min=rmxl_min)

    z_idx = torch.arange(nz + 1, device=mxl.device, dtype=mxl.dtype) * dz
    x = z_idx.view(1, nz + 1, 1)

    # --- forward sweep (surface -> bottom) ---
    mxl_fwd = mxl.clone()
    mxl_fwd[:, 0, :] = cfg.tke_rmxl0 + cfg.Hmin
    lup = torch.cummin(mxl_fwd - x, dim=1)[0] + x

    # --- backward sweep (bottom -> surface) ---
    mxl_bwd = mxl.clone()
    mxl_bwd[:, -1, :] = rmxl_min                # bottom BC
    mxl_bwd = torch.flip(mxl_bwd, dims=[1])
    x_flip  = torch.flip(x, dims=[1])
    ldn = torch.flip(torch.cummin(mxl_bwd + x_flip, dim=1)[0], dims=[1]) - x

    le = torch.sqrt(lup * ldn).clamp(min=rmxl_min)
    lk = torch.minimum(lup, ldn).clamp(min=rmxl_min)
    return lk, le