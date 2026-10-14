import torch
from taia1dvfd.core import Model, pcr_solve

def temperature_equation_init(model: Model, temp_init: torch.Tensor | None = None):
    zc = model.grid.z_c
    cfg = model.config

    # --- CAREFULL, FORCING LOTUS ONLY --- 
    if temp_init is not None:
        T0 = temp_init.clone()
        
    else:
        T0 = torch.zeros_like(model.state.temp)
        dT = cfg.Tsurf - cfg.Tbotm
        T0[:] = cfg.Tsurf - dT * (-zc - cfg.Hmin) / (cfg.Hmax - cfg.Hmin)
        
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
    sig_z    = state.sig_z_t
    
    # --- Parameters ---
    fkw = cfg.farr_kwargs
    dz = grd.dz
    dt = cfg.dt
    nz = cfg.nz
    
    # --- Diffusion operator ---
    d_aztdz = ops.d_aztdz.clone()

    # --- Neumann BC ---
    if cfg.temp_bc_top == "neumann":
        azt_dn = azt[:, 1, 0] # face below surface cell
        d_aztdz[:, 0, 1] =  azt_dn / dz**2
        d_aztdz[:, 0, 0] = -azt_dn / dz**2
    
    if cfg.temp_bc_btm == "neumann":
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
    rhs = temp_old + Fk*dt
    
    # if cfg.temp_stochastic:
    #     Gk = sig_z * (Dz @ temp_old)
    #     rhs -= Gk*state.dB_c
    
    # --- Stochastic transport (face-based, discrete fluctuation-dissipation balance) ---
    if cfg.temp_stochastic:
        dTdz_w  = ops.d_dz_c2w @ temp_old                  # (ne, nz+1, 1) dT/dz at faces, 0 at boundary faces
        noise_w = state.sig_z_t * dTdz_w * state.dB_w      # sigma_t and dB at faces
        Gk      = (noise_w[:, :-1] + noise_w[:, 1:]) / 2**0.5
        rhs    -= Gk
    
    if cfg.temp_bc_top == "neumann":
        rhs[:,  0, 0] += -dt * Fsurf / dz
    elif cfg.temp_bc_top == "dirichlet":  
            Lop[:, 0, :] = 0.0
            Lop[:, 0, 0] = 1.0
            rhs[:, 0, :] = temp_surf.unsqueeze(-1)
    
    if cfg.temp_bc_btm == "neumann":
        rhs[:, -1, 0] += dt * Fbotm / dz
    elif cfg.temp_bc_btm == "dirichlet":
        Lop[:, -1,  :] = 0.0
        Lop[:, -1, -1] = 1.0
        rhs[:, -1,  :] = temp_botm.unsqueeze(-1)
    
    # --- Solver ---
    # state.temp[:] = torch.linalg.solve(Lop, rhs)
    # state.temp[:] = thomas_solve(Lop, rhs)
    state.temp[:] = pcr_solve(Lop, rhs)