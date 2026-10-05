import torch
from taia1dvfd.core import Model

def momentum_equation_init(model: Model) -> None:
    pass

def momentum_equation_step(model: Model) -> None:
    state = model.state
    cfg   = model.config
    grd   = model.grid
    ops   = model.operators
    
    # --- Variables ---
    ue = state.ue
    us = state.us
    ua = state.ua
    azm   = state.az_m
    sig_z = state.sig_z
    sig_h = state.sig_h
    
    # --- Parameters ---
    fkw = cfg.farr_kwargs
    ckw = cfg.carr_kwargs
    dz = grd.dz
    dt = cfg.dt
    sqrt_dt = cfg.sqrt_dt
    cdrho = cfg.cdrho
    nz = cfg.nz
    f = cfg.f
    Dz = ops.d_dz_w
    Idc = torch.eye(nz+1, **ckw)

    # --- Neumann BC ---
    d_azmdz = ops.d_azmdz.clone()
    azm_dn = azm[:, 1, 0]
    d_azmdz[:, 0, 1] =  azm_dn / dz**2
    d_azmdz[:, 0, 0] = -azm_dn / dz**2

    azm_up = azm[:, nz-1, 0]
    d_azmdz[:, -1, -2] =  azm_up / dz**2
    d_azmdz[:, -1, -1] = -azm_up / dz**2

    if cfg.wave_mixing:
        # --- Linear operator ---
        A = d_azmdz - 1j * f * Idc
        L = Idc - dt * A

        # --- RHS --- #
        Fk = A @ us
        rhs = ue + Fk * dt

        if cfg.vel_stochastic:
            dB = sqrt_dt * torch.randn_like(ue, **fkw)
            Gk = 1j * f * sig_h + (sig_z * Dz) @ (ue + us)
            rhs -= Gk*dB
    else:
        # --- Linear operator ---
        L = Idc - dt * (d_azmdz - 1j*f*Idc)
        
        # --- RHS --- #
        Fk = -1j*f*us
        rhs = ue + Fk * dt
        
        if cfg.vel_stochastic:
            dB = sqrt_dt * torch.randn_like(ue, **fkw)
            Gk = 1j * f * sig_h + (sig_z * Dz) @ ue
            rhs -= Gk*dB

    # --- Boundary forcing ---
    tau_w = (cdrho * torch.abs(ua) * ua).squeeze(-1)
    tau_s = (azm[:, 0, :] * state.us[:, 0, :] / cfg.hs)
    
    Fsurf = tau_w + tau_s if cfg.wave_mixing else tau_w
    Fbotm = torch.zeros_like(tau_w)
    
    rhs[:,  0, :] += dt * Fsurf / dz
    rhs[:, -1, :] += dt * Fbotm / dz

    # Solve L ue^{n+1} = RHS
    ue_new = torch.linalg.solve(L, rhs)

    model.state.ue[:] = ue_new