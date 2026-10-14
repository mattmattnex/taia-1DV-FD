import torch

from taia1dvfd.core  import Config, Grid, State, Model, Operators
from taia1dvfd.utils import SaveNetCDF, StepProfiler

import taia1dvfd.physics as phy
import taia1dvfd.mixing  as mix
import taia1dvfd.equation_of_state as eos


for model_conf in ['PDE', 'SPDE-noWM', 'SPDE']:
    cfg = Config()

    cfg.ne   = 1000
    cfg.nz   = 300
    cfg.Hmin = 0.0
    cfg.Hmax = 100.0
    cfg.lat  = 34.98
    cfg.device = "cuda:0"

    # --- Time ---
    cfg.dt   = 1800
    cfg.tend = 3600.0 * 24 * 30

    # --- Wind ---
    cfg.wind_mean = 6.0 + 0.0j
    cfg.wind_std  = 0.0 + 1.0j

    # --- Config bool ---
    if model_conf == 'PDE':
        cfg.wave_mixing = False
        cfg.vel_stochastic  = False
    
    elif model_conf == 'SPDE-noWM':
        cfg.wave_mixing = False
        cfg.vel_stochastic  = True
    
    elif model_conf == 'SPDE':
        cfg.wave_mixing = True
        cfg.vel_stochastic  = True

    # ---------------------------------
    # ----- MIXING model choice -------
    # ---------------------------------
    cfg.kpp_regularization = False
    cfg.kpp_azm_min = 1e-4
    cfg.kpp_azt_min = 1e-4
                
    # --- Finalization ---
    cfg.save_folder = "/home/mnex/Desktop/1DV_FD/experiments/0_tests/output/"
    cfg.file_name = f"longtest_{model_conf}.nc"
    cfg.var_to_save = ['ue']
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
    mix.kpp_init(model)

    savefile.init_file(model.state)
    savefile.save(model.state)
            
    # ================== Main Loop ==================
    # Order : Wind -> Mixing + Ops -> Noise -> Momentum
    # ===============================================
    for kt in range(cfg.nt):
        prof.start()
        
        # --- Wind forcing ---
        phy.wind_ornstein_step(model)
        u_star = phy.compute_ustar(model).squeeze(-1)
        
        # --- Vertical Mixing ---
        mix.kpp_step(u_star, model)
        operators.update_momentum_diffusion(state.az_m)
        operators.update_tracers_diffusion(state.az_t)
        prof.mark("Mix")

        # --- Noise ---
        phy.noise_LU_step(model)
        prof.mark("Sig")
        
        # --- Momentum --- 
        phy.momentum_equation_step(model)
        prof.mark("Phy")

        # --- Update time ---
        model.state.time += model.config.dt
        
        # --- Save file ---
        savefile.save(model.state)

        if kt%5 == 0:
            print(prof.format_step(kt, cfg.nt))

    print(prof.format_summary())
    savefile.close()