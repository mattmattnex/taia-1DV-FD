import torch
from taia1dvfd.core import Model

def compute_ustar(model: Model):
    cfg = model.config
    ua = model.state.ua

    return cfg.sqrt_cdrho * torch.abs(ua)

def wind_ornstein_init(model: Model):
    model.state.ua[:] = model.config.wind_mean

def wind_ornstein_step(model: Model):
    cfg = model.config
    ckw = cfg.carr_kwargs
    
    # Fetch variables
    ua = model.state.ua
    
    # Fetch params
    ua_m = cfg.wind_mean
    
    # Wind forcing coefficients
    c1 = torch.exp(torch.tensor(-cfg.dt / cfg.wind_Ta, **cfg.farr_kwargs))
    c2 = torch.sqrt(1 - c1**2) * cfg.wind_std
    
    # Random wind fluctuation
    rd_cplx = torch.randn_like(ua, **ckw)

    # Wind velocity update
    ua_new = ua_m + c1 * (ua - ua_m) + c2 * rd_cplx
    
    model.state.ua[:] = ua_new