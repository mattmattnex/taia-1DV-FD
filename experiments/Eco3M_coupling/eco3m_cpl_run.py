import torch
import numpy as np

from taia1dvfd.core  import Config, Grid, State, Model, Operators
from taia1dvfd.utils import SaveNetCDF, StepProfiler

import taia1dvfd.physics as phy
import taia1dvfd.mixing  as mix
import taia1dvfd.equation_of_state as eos
import taia1dvfd.biogeochem as bio

cfg = Config()

cfg.ne   = 100
cfg.nz   = 256
cfg.Hmin = 1.0
cfg.Hmax = 100.0
cfg.lat  = 34.98
cfg.device = "cuda:0"

# --- Time ---
cfg.dt   = 1800
cfg.tend = 3600.0 * 24 * 30

# --- Wind ---
cfg.wind_mean = 5.0 + 0.0j
cfg.wind_std  = 1.0

# --- Wave ---
cfg.wave_mixing = True
cfg.tke_wave_mix = False
cfg.gls_wave_mix = False

# --- Stochastic bool ---
cfg.temp_stochastic = False
cfg.salt_stochastic = False
cfg.vel_stochastic  = True

# --- EOS ---
cfg.Tsurf  = 20
cfg.Tbotm  = 20
cfg.Ssurf = 38
cfg.Sbotm = 38

# --- TKE model ---
cfg.tke_pdl = True        
cfg.tke_ebb = 67.83
cfg.tke_azm_min = 1.2e-4
cfg.tke_azt_min = 1.2e-5

# --- Finalization ---
cfg.save_folder = "/home/mnex/Desktop/1DV_FD/experiments/Eco3M_coupling/output/"
cfg.file_name = f"eco3M_ne{cfg.ne}.nc"
cfg.var_to_save = ['ue', 'mld', 'N2', 'temp']
cfg.finalize()

# --- Model objects ---
grid   = Grid(cfg)
state  = State(cfg)
operators = Operators(cfg, grid)
model    = Model(cfg, grid, state, operators)
savefile = SaveNetCDF(model)
prof    = StepProfiler()

# --- Initializations ---
phy.noise_LU_init(model)
phy.wind_ornstein_init(model)

mix.tke_init(model)
mix.mld_init(model)

eos.temperature_equation_init(model)
eos.salinity_equation_init(model)
eos.eos_init(model)

# --- Initialize Eco3M ---
Eco3M = bio.bgc_coupling_init(model)
Eco3M.CONC[:, 0] = 1.0e-7
Eco3M.CONC[:, 1] = 5.0e-6

savefile.init_file(model.state)
savefile.save(model.state)

# ================== Main Loop ==================
for kt in range(cfg.nt):
    prof.start()
    
    # --- Wind forcing ---
    phy.wind_ornstein_step(model)
    u_star = phy.compute_ustar(model).squeeze(-1)
    prof.mark("Wind")
    
    # --- Vertical Mixing ---
    mix.tke_step(u_star, model)
    prof.mark("Mix")
    
        # --- Update mixing operators ---
    operators.update_momentum_diffusion(state.az_m)
    operators.update_tracers_diffusion(state.az_t)
    prof.mark("Ops")
    
    # --- EOS ---
    eos.temperature_equation_step(model)
    eos.salinity_equation_step(model)
    eos.eos_step(model)
    prof.mark("eos")

    # --- MLD ---
    mix.mld_step(model)
    prof.mark("mld")

    # --- Physics and stochastic
    phy.noise_LU_step(model)        
    phy.momentum_equation_step(model)
    prof.mark("Phy")
    
    # --- Eco3M ---
    bio.bgc_coupling_step(model, Eco3M)
    
    # --- Update time ---
    model.state.time += model.config.dt
    
    # --- Save file ---
    savefile.save(model.state)

    if kt%5 == 0:
        print(prof.format_step(kt, cfg.nt))

print(prof.format_summary())
savefile.close()