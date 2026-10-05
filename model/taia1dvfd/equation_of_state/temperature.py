import torch
from taia1dvfd.core import Model

def temperature_equation_init(model: Model, temp_init: torch.Tensor | None = None):
    zc = model.grid.z_c
    cfg = model.config

    # --- CAREFULL, FORCING LOTUS ONLY --- 
    if temp_init is not None:
        T0 = temp_init.clone()
        
    else:
        T0 = torch.zeros_like(model.state.temp)
        dT = cfg.Tsurf - cfg.Tbotm
        T0[:] = cfg.Tsurf + dT * zc / model.config.Hmax
        
    model.state.temp[:] = T0
    
def temperature_equation_step(model: Model, 
                              temp_surf: torch.Tensor | None = None,
                              temp_botm: torch.Tensor | None = None,
                              ):
    state = model.state
    cfg   = model.config
    grd   = model.grid
    ops   = model.operators
    
    # --- Variables ---
    temp_old = state.temp
    azt      = state.az_t
    sig_z    = state.sig_z
    
    # --- Parameters ---
    fkw = cfg.farr_kwargs
    sqrt_dt = cfg.sqrt_dt
    dz = grd.dz
    dt = cfg.dt
    nz = cfg.nz
    
    bc_top = "dirichlet"
    bc_btm = "neumann"
    Fsurf = 0.0
    Fbotm = 0.0
    
    # --- Diffusion operator ---
    d_aztdz = ops.d_aztdz.clone()

    # --- Neumann BC ---
    if bc_top == "neumann":
        azt_dn = azt[:, 1, 0] # face below surface cell
        d_aztdz[:, 0, 1] =  azt_dn / dz**2
        d_aztdz[:, 0, 0] = -azt_dn / dz**2
    
    if bc_btm == "neumann":
        azt_up = azt[:, nz-1, 0] # face above bottom cell
        d_aztdz[:, -1, -2] =  azt_up / dz**2
        d_aztdz[:, -1, -1] = -azt_up / dz**2

    # --- Linear operator ----
    Id = torch.eye(nz, **fkw)
    Lop = Id - dt*d_aztdz
    Dz = ops.d_dz_c
        
    # --- RHS ---
    Fk = 0
    rhs = temp_old + Fk*dt
    
    if cfg.temp_stochastic:
        dB = sqrt_dt * torch.randn_like(temp_old, **fkw)
        Gk = grd.to_c(sig_z) * (Dz @ temp_old)
        rhs -= Gk*dB
    
    if bc_top == "neumann":
        rhs[:,  0, 0] += -dt * Fsurf / dz
    
    if bc_btm == "neumann":
        rhs[:, -1, 0] += -dt * Fbotm / dz

    elif bc_top == "dirichlet":  
        Lop[:, 0, :] = 0.0
        Lop[:, 0, 0] = 1.0
        rhs[:, 0, :] = temp_surf.unsqueeze(-1)
    
    elif bc_btm == "dirichlet":
        Lop[:, -1,  :] = 0.0
        Lop[:, -1, -1] = 1.0
        rhs[:, -1,  :] = temp_botm.unsqueeze(-1)
    
    # --- Solver ---
    state.temp[:] = torch.linalg.solve(Lop, rhs)