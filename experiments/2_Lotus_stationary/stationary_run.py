import torch

from taia1dvfd.core  import Config, Grid, State, Model, Operators
from taia1dvfd.utils import SaveNetCDF, StepProfiler

import taia1dvfd.physics as phy
import taia1dvfd.mixing  as mix
import taia1dvfd.equation_of_state as eos

CLOSURES = [
    ('kpp', 'kpp'),
    ('gls', 'k-kl'),
    ('gls', 'k-epsilon'),
    ('gls', 'k-omega'),
    ('gls', 'gen'),
    ('tke', 'tke')
]

for prefix, mix_model in CLOSURES:
    cfg = Config()

    cfg.ne   = 1000
    cfg.nz   = 300
    cfg.Hmin = 1.0
    cfg.Hmax = 100.0
    cfg.lat  = 34.98
    cfg.device = "cuda:0"

    # --- Time ---
    cfg.dt   = 3600
    cfg.tend = 3600.0 * 24 * 30

    # --- Wind ---
    cfg.wind_mean = 6.0 + 0.0j
    cfg.wind_std  = 0.0 + 1.0j

    # --- Wave ---
    cfg.wave_mixing = True
    
    # --- Stochastic bool ---
    cfg.temp_stochastic = False
    cfg.salt_stochastic = False
    cfg.vel_stochastic  = True

    # --- EOS ---
    cfg.Tsurf  = 20
    cfg.Tbotm  = 20
    cfg.Ssurf = 38
    cfg.Sbotm = 38
    
    # -------------------------
    # ----- MIXING model ------
    # -------------------------
    # Params for KPP
    cfg.kpp_regularization = False
    cfg.kpp_azm_min = 1e-4
    cfg.kpp_azt_min = 1e-4
    
    # Params for tke closure
    cfg.tke_wave_mix = True
    cfg.tke_pdl = True        
    cfg.tke_ebb = 67.83
    cfg.tke_azm_min = 1e-4
    cfg.tke_azt_min = 1e-4
    
    # Params for GLS
    cfg.gls_stabf = 'CA'
    cfg.gls_model = mix_model
    cfg.gls_azm_min  = 1e-4
    cfg.gls_azt_min  = 1e-4
    cfg.gls_tke_min  = 1e-8
    cfg.gls_diss_min = 1e-12
    
    # --- Finalization ---
    cfg.save_folder = "/home/mnex/Desktop/1DV_FD/experiments/LOTUS/output/"
    cfg.file_name =  f'lotus_closure_{mix_model}_ne{cfg.ne}.nc'
    cfg.var_to_save = ['ue', 'az_m', 'az_t', 'temp']
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

    if mix_model == 'kpp':
        mix.kpp_init(model)
    elif mix_model == 'tke':
        mix.tke_init(model)
    elif mix_model in ['k-kl', 'k-epsilon', 'k-omega', 'gen']:
        mix.gls_init(model)
    else:
        raise ValueError('Mix model not correct')
        
    mix.mld_init(model)

    eos.temperature_equation_init(model)
    eos.salinity_equation_init(model)
    eos.eos_init(model)

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
        if mix_model == 'kpp':
            mix.kpp_step(u_star, model)
        elif mix_model == 'tke':
            mix.tke_step(u_star, model)
        else:
            mix.gls_step(u_star, model)
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
        
        # --- Update time ---
        model.state.time += model.config.dt
        
        # --- Save file ---
        savefile.save(model.state)

        if kt%5 == 0:
            print(prof.format_step(kt, cfg.nt))

    print(prof.format_summary())
    savefile.close()