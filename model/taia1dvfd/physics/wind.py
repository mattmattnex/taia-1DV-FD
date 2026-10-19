import torch
import math
from taia1dvfd.core import Model

def compute_ustar(model: Model):
    cfg = model.config
    ua = model.state.ua

    return cfg.sqrt_cdrho * torch.abs(ua)

def wind_ornstein_init(model: Model):
    model.state.ua[:] = model.config.wind_mean

def wind_ornstein_step(model: Model):
    cfg = model.config
    
    # Fetch variables
    ua = model.state.ua
    
    # Fetch params
    ua_m = cfg.wind_mean
    
    # Wind forcing coefficients
    c1 = c1 = torch.exp(-torch.tensor(cfg.dt / cfg.wind_Ta, **cfg.farr_kwargs))
    c2 = torch.sqrt(1 - c1**2) * torch.tensor(cfg.wind_std, **cfg.farr_kwargs)
    
    # Random wind fluctuation
    rn = torch.randn(2, cfg.ne, 1, 1, **cfg.farr_kwargs)

    # Wind velocity update
    ua_new = ua_m + c1 * (ua - ua_m) + c2[0]*rn[0] + 1j*c2[1]*rn[1]
    
    model.state.ua[:] = ua_new