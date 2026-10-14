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
    cfg   = model.config
    state = model.state
    Dcw   = model.operators.d_dz_c2w

    # Buoyancy anomaly b = -g (rho - rho_o) / rho_o  (no float32 cancellation around 1000)
    buoy = cfg.g * (cfg.beta_T * (state.temp - cfg.eos_T0)
                  - cfg.beta_S * (state.salt - cfg.eos_S0))

    N2 = Dcw @ buoy                    # N2 = db/dz with the true d/dz operator
    N2[:, 0,  :] = N2[:,  1, :]
    N2[:, -1, :] = N2[:, -2, :]

    state.N2[:] = N2