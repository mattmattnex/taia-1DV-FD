import torch
from taia1dvfd.core import Model, pcr_solve

def salinity_equation_init(model: Model, salt_init: torch.Tensor | None = None):
    zc = model.grid.z_c
    cfg = model.config

    # --- CAREFULL, FORCING LOTUS ONLY --- 
    if salt_init is not None:
        S0 = salt_init.clone()
        
    else:
        S0 = torch.zeros_like(model.state.salt)
        dS = cfg.Ssurf - cfg.Sbotm
        S0[:] = cfg.Ssurf - dS * (-zc - cfg.Hmin) / (cfg.Hmax - cfg.Hmin)
        
    model.state.salt[:] = S0

def salinity_equation_step(model: Model, 
                              salt_surf: torch.Tensor | None = None,
                              salt_botm: torch.Tensor | None = None,
                              ):
    state = model.state
    cfg   = model.config
    grd   = model.grid
    ops   = model.operators
    
    # --- Variables ---
    salt_old = state.salt
    azt      = state.az_t
    sig_z    = state.sig_z_t
    
    # --- Parameters ---
    fkw = cfg.farr_kwargs
    dz = grd.dz
    dt = cfg.dt
    nz = cfg.nz
        
    # --- Diffusion operator ---
    d_aztdz = ops.d_aztdz.clone()

    # --- Neumann BC ---
    if cfg.salt_bc_top == "neumann":
        azt_dn = azt[:, 1, 0] # face below surface cell
        d_aztdz[:, 0, 1] =  azt_dn / dz**2
        d_aztdz[:, 0, 0] = -azt_dn / dz**2
    
    if cfg.salt_bc_btm == "neumann":
        azt_up = azt[:, nz-1, 0] # face above bottom cell
        d_aztdz[:, -1, -2] =  azt_up / dz**2
        d_aztdz[:, -1, -1] = -azt_up / dz**2

    # --- Linear operator ----
    Id = torch.eye(nz, **fkw)
    Lop = Id - dt*d_aztdz
    Dz = ops.d_dz_c

    # --- RHS ---
    Fsurf = 0.0
    Fbotm = 0.0
    Fk = 0
    rhs = salt_old + Fk*dt
    
    if cfg.salt_stochastic:
        Gk = sig_z * (Dz @ salt_old)
        rhs -= Gk * state.dB_c
    
    if cfg.salt_bc_top == "neumann":
        rhs[:,  0, 0] += -dt * Fsurf / dz
    elif cfg.salt_bc_top == "dirichlet":  
            Lop[:, 0, :] = 0.0
            Lop[:, 0, 0] = 1.0
            rhs[:, 0, :] = salt_surf.unsqueeze(-1)
    
    if cfg.salt_bc_btm == "neumann":
        rhs[:, -1, 0] += dt * Fbotm / dz
    elif cfg.salt_bc_btm == "dirichlet":
        Lop[:, -1,  :] = 0.0
        Lop[:, -1, -1] = 1.0
        rhs[:, -1,  :] = salt_botm.unsqueeze(-1)
    
    # --- Solver ---
    # state.salt[:] = torch.linalg.solve(Lop, rhs)
    # state.salt[:] = thomas_solve(Lop, rhs)
    state.salt[:] = pcr_solve(Lop, rhs)