import torch
from taia1dvfd.core import Model
from .tracer import tracer_equation_step

def temperature_equation_init(model: Model, temp_init: torch.Tensor | None = None):
    zc = model.grid.z_c
    cfg = model.config

    # --- CAREFULL, FORCING LOTUS ONLY ---
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
    cfg = model.config

    model.state.temp[:] = tracer_equation_step(
        model,
        tracer     = model.state.temp,
        tracer_ref = cfg.eos_T0,
        stochastic = cfg.temp_stochastic,
        bc_top     = cfg.temp_bc_top,
        bc_btm     = cfg.temp_bc_btm,
        value_top  = temp_surf,
        value_btm  = temp_botm,
        flux_top   = 0.0,          # w'T'(0), positive upward (heat flux: -Q / (rho_o cp))
        flux_btm   = 0.0,
    )
