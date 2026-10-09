import torch
from taia1dvfd.core import Model

def eos_init(model: Model):
    eos_step(model)

def eos_step(model: Model):
    _density(model)
    _brunt_vaisala(model)
    
def _density(model: Model):
    cfg = model.config
    state = model.state

    T = state.temp
    S = state.salt
    
    model.state.rho[:] = cfg.rho_o * ((1 - cfg.beta_T * (T - cfg.eos_T0)) \
                                         + cfg.beta_S * (S - cfg.eos_S0))
    
def _brunt_vaisala(model: Model):
    cfg = model.config
    ops = model.operators

    rho = model.state.rho
    Dcw = ops.d_dz_c2w

    # Derivative center to faces nz -> nz+1
    N2 = -cfg.g / cfg.rho_o * (Dcw @ rho)
    
    # zero-gradient BC
    N2[:, 0,  :] = N2[:,  1, :]
    N2[:, -1, :] = N2[:, -2, :]

    model.state.N2[:]  = N2