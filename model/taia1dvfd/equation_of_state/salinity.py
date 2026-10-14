import torch
from taia1dvfd.core import Model
from .tracer import tracer_equation_step

def salinity_equation_init(model: Model, salt_init: torch.Tensor | None = None):
    zc = model.grid.z_c
    cfg = model.config

    # --- CAREFULL, FORCING LOTUS ONLY ---
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
    cfg = model.config

    model.state.salt[:] = tracer_equation_step(
        model,
        tracer     = model.state.salt,
        tracer_ref = cfg.eos_S0,
        stochastic = cfg.salt_stochastic,
        bc_top     = cfg.salt_bc_top,
        bc_btm     = cfg.salt_bc_btm,
        value_top  = salt_surf,
        value_btm  = salt_botm,
        flux_top   = 0.0,          # w'S'(0), positive upward (freshwater: -S0 (E - P))
        flux_btm   = 0.0,
    )
