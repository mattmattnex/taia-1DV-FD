import torch
from taia1dvfd.core import Model

def noise_LU_init(model: Model):
    cfg  = model.config
    grid = model.grid
    
    # Wave amplitude and direction
    us_amp = cfg.us0 * torch.exp(grid.z_w / cfg.hs + 1j * cfg.theta_mean)
    angle_rd = cfg.theta_std * torch.randn(cfg.ne, 1, 1, **cfg.farr_kwargs)
    us = us_amp * torch.exp(1j * angle_rd)
    
    # Set wave 
    model.state.us[:] = us

def noise_LU_step(model:Model):
    cfg = model.config
    state = model.state
    
    # Fetch
    us   = model.state.us
    az_m = model.state.az_m
        
    # Compute noise
    sz = torch.sqrt(2 * az_m)
    sh = 2 * cfg.hs * us  / sz

    # Update
    state.sig_z[:] = sz
    state.sig_h[:] = sh
    state.dB_w[:] = cfg.sqrt_dt * torch.randn(cfg.ne, cfg.nz+1, 1, **cfg.farr_kwargs)
    state.dB_c[:] = (state.dB_w[:, :-1] + state.dB_w[:, 1:]) / 2**0.5 # unit variance, correlated with both faces, design choice