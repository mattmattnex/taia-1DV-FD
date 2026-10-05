import torch
from taia1dvfd.core import Model

def kpp_init(model: Model):
    cfg = model.config
    state = model.state
    
    state.az_m[:] = cfg.kpp_azm_min
    state.az_t[:] = cfg.kpp_azt_min

def kpp_step(u_star, model:Model):
    cfg = model.config
    
    # Fetch params
    z  = model.grid.z_w
    
    # Coefficients for kpp calculations
    a1 = cfg.kpp_c2 / cfg.f

    # Boundary layer depth
    u_star = u_star.unsqueeze(-1)
    h = a1 * u_star
    blMask = torch.where(abs(z) <= h, 1., 0.)
    zt = - z / h

    # Non-dimensional shear
    G_zt = zt * (1 - zt)**2
    
    if cfg.kpp_regularization:
        zt_0 = cfg.kpp_zt_0
        G_zt += torch.heaviside(zt_0 - zt, zt[-1]) * (zt - zt_0)**2 / (2 * zt_0)
    
    #  Eddy viscosity
    a2 = cfg.kpp_c1 * u_star * h * G_zt
    nu = blMask * a2
        
    model.state.az_m[:] = nu + cfg.kpp_azm_min
    model.state.az_t[:] = nu + cfg.kpp_azt_min