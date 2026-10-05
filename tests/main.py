import time
import torch

from taia1dvfd.core  import Config, Grid, State, Model, Operators
from taia1dvfd.utils import Animator, StepProfiler

import taia1dvfd.physics as phy
import taia1dvfd.mixing  as mix
import taia1dvfd.equation_of_state as eos

cfg = Config()

cfg.ne   = 20
cfg.nz   = 256
cfg.Hmin = 1.0
cfg.Hmax = 50.0
cfg.lat  = 0.0
cfg.device = "cuda:0"

# --- Time ---
cfg.dt   = 36
cfg.tend = 3600.0 * 30

# --- Stochastic bool ---
cfg.temp_stochastic = True
cfg.salt_stochastic = True
cfg.vel_stochastic  = True

# --- EOS ---
cfg.Tsurf  = 13
cfg.Tbotm  = cfg.Tsurf - (1e-4 / (cfg.g * cfg.beta_T)) * (cfg.Hmax - cfg.Hmin)
cfg.Ssurf = 38
cfg.Sbotm = 38

# --- TKE model ---
cfg.azm_min = 1.2e-4
cfg.azt_min = 1.2e-4
cfg.tke_ebb = 67.83 # 3.75 else no wave
cfg.tke_pdl  = True

# --- Finalization ---
# cfg.save_folder = "/home/mnex/Desktop/1DV_model_PHD/experiments/Kato_Phillips/output/tke/"
# cfg.file_name = f"K-P_tke_{cfg.dt}_{cfg.nz}.nc"
# cfg.var_to_save = ['temp', 'mld', 'N2']
cfg.finalize()


# --- Model objects ---
grid   = Grid(cfg)
state  = State(cfg)
operators = Operators(cfg, grid)
model    = Model(cfg, grid, state, operators)
anim    = Animator(model)
prof    = StepProfiler()

# --- Init functions ---
phy.wind_ornstein_init(model)
phy.noise_LU_init(model)

# mix.kpp_init(model)
mix.tke_init(model)

eos.eos_init(model)
eos.temperature_equation_init(model)
eos.salinity_equation_init(model)
phy.momentum_equation_init(model)

# ===============================================
# ================== MAIN LOOP ==================

for kt in range(cfg.nt):
    prof.start()

    # --- Wind ---
    phy.wind_ornstein_step(model)
    u_star = phy.compute_ustar(model).squeeze(-1)
    prof.mark("Wind")
    
    # --- Mixing ---
    # mix.kpp_step(u_star, model)
    mix.tke_step(u_star, model)
    prof.mark("Mixing")
    
    # --- Update mixing operators ---
    operators.update_momentum_diffusion(state.az_m)
    operators.update_tracers_diffusion(state.az_t)
    prof.mark("Ops")
    
    # --- EOS ---
    eos.temperature_equation_step(model)
    eos.salinity_equation_step(model)
    eos.eos_step(model)
    prof.mark("eos")
    
    # --- Physics ---
    phy.noise_LU_step(model)
    phy.momentum_equation_step(model)
    prof.mark("Phy")
    
    # --- Update time ---
    model.state.time += model.config.dt

    # --- Animation ---
    if kt%5 == 0:
        print(prof.format_step(kt, cfg.nt))
        anim.update_plot(model)