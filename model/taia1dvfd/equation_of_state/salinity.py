import torch
from taia1dvfd.core import Model

def salinity_equation_init(model: Model):
    z = model.grid.z_c
    cfg = model.config
    
    dS = cfg.Ssurf - cfg.Sbotm
    salt_init = cfg.Ssurf + dS * z / model.config.Hmax

    model.state.salt[:] = salt_init

def salinity_equation_step(model: Model):
    state = model.state
    cfg   = model.config
    grd   = model.grid
    ops   = model.operators
    
    # --- Variables ---
    salt_old = state.salt
    azt      = state.az_t
    sig_z    = state.sig_z
    
    # --- Parameters ---
    fkw = cfg.farr_kwargs
    sqrt_dt = cfg.sqrt_dt
    dz = grd.dz
    dt = cfg.dt
    nz = cfg.nz
    
    bc = 'neumann'
    Fsurf = 0.0
    Fbotm = 0.0
    
    # --- Diffusion operator ---
    d_aztdz = ops.d_aztdz.clone()

    # --- Neumann BC ---
    if bc == "neumann":
        azt_dn = azt[:, 1, 0] # face below surface cell
        d_aztdz[:, 0, 1] =  azt_dn / dz**2
        d_aztdz[:, 0, 0] = -azt_dn / dz**2

        azt_up = azt[:, nz-1, 0] # face above bottom cell
        d_aztdz[:, -1, -2] =  azt_up / dz**2
        d_aztdz[:, -1, -1] = -azt_up / dz**2

    # --- Linear operator ----
    Id = torch.eye(nz, **fkw)
    Lop = Id - dt*d_aztdz
    Dz = ops.d_dz_c
        
    # --- RHS ---
    Fk = 0
    rhs = salt_old + Fk*dt
    
    if cfg.salt_stochastic:
        dB = sqrt_dt * torch.randn_like(salt_old, **fkw)
        Gk = grd.to_c(sig_z) * (Dz @ salt_old)
        rhs -= Gk*dB
    
    if bc == "neumann":
        rhs[:,  0, 0] += -dt * Fsurf / dz
        rhs[:, -1, 0] += -dt * Fbotm / dz

    elif bc == "dirichlet":  
        Lop[:, 0, :] = 0.0
        Lop[:, 0, 0] = 1.0
        rhs[:, 0, :] = 0.5 # Surface value
        
        Lop[:, -1,  :] = 0.0
        Lop[:, -1, -1] = 1.0
        rhs[:, -1,  :] = 0.5 # Bottom value    
    
    # --- Solver ---
    state.salt[:] = torch.linalg.solve(Lop, rhs)