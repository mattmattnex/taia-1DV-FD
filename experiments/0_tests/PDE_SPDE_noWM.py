import torch

from taia1dvfd.core  import Config, Grid, State, Model
from taia1dvfd.utils import SaveNetCDF, StepProfiler

import taia1dvfd.physics as phy
import taia1dvfd.mixing  as mix
import taia1dvfd.equation_of_state as eos

SEED = 1234
CONFIGS = ['PDE', 'SPDE-noWM', 'SPDE']
NE   = 1000
NZ   = 300
HMIN = 1.0
HMAX = 100.0
LAT  = 34.98
DT   = 3600.0
TEND = 3600.0 * 24 * 30

for conf in CONFIGS:
    print(f'Running : {conf}')

    # --- Same random stream for every configuration ---
    torch.manual_seed(SEED)

    cfg = Config()
    cfg.ne   = NE
    cfg.nz   = NZ
    cfg.Hmin = HMIN
    cfg.Hmax = HMAX
    cfg.lat  = LAT
    cfg.device = "cuda:0"

    # --- Time ---
    cfg.dt   = DT
    cfg.tend = TEND

    # --- Wind ---
    cfg.wind_mean = 6.0 + 0.0j
    cfg.wind_std  = 0.0 + 1.0j
    
    # --- Stochastic bool ---
    cfg.temp_stochastic = False
    cfg.salt_stochastic = False

    if conf == 'PDE':
        cfg.vel_stochastic = False
        cfg.wave_mixing    = False
    elif conf == 'SPDE-noWM':
            cfg.vel_stochastic = True
            cfg.wave_mixing    = False
    elif conf == 'SPDE':
            cfg.vel_stochastic = True
            cfg.wave_mixing    = True
    else: 
        raise ValueError("Config not valid")

    # --- EOS: neutral column (uniform T and S, N2 = 0) ---
    cfg.Tsurf = cfg.Tbotm = 20.0
    cfg.Ssurf = cfg.Sbotm = 38.0
    cfg.temp_bc_top = cfg.temp_bc_btm = "neumann"
    cfg.salt_bc_top = cfg.salt_bc_btm = "neumann"
    
    # ---------------------------------
    # ----- MIXING model choice -------
    # ---------------------------------
    cfg.kpp_regularization = False
    cfg.kpp_azm_min = 1e-4
    cfg.kpp_azt_min = 1e-4

# --- Finalization ---
    cfg.save_folder = f"/home/mnex/Desktop/1DV_FD/experiments/0_tests/output/"
    cfg.file_name   = f"{conf}_ne{cfg.ne}.nc"
    cfg.var_to_save = ['ue', 'ua', 'az_m', 'tke']
    cfg.finalize()

    # --- Model objects ---
    grid     = Grid(cfg)
    state    = State(cfg)
    model    = Model(cfg, grid, state)
    savefile = SaveNetCDF(model)
    prof     = StepProfiler()

    # --- Initializations ---
    phy.noise_LU_init(model)            # Stokes drift with random direction
    phy.wind_ornstein_init(model)
    
    mix.kpp_init(model)
    
    eos.temperature_equation_init(model)
    eos.salinity_equation_init(model)
    eos.eos_init(model)                 # N2 = 0 exactly (uniform anomalies)

    # --- Output file ---
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
        prof.mark("Mix")

        # --- Noise ---
        phy.noise_LU_step(model)
        prof.mark("Sig")

        # --- Temp Sal ---
        eos.temperature_equation_step(model)
        eos.salinity_equation_step(model)
        prof.mark("T-S")

        # --- Momentum ---
        phy.momentum_equation_step(model)
        prof.mark("Phy")

        # --- EOS ---
        eos.eos_step(model)

        # --- Update time ---
        model.state.time += model.config.dt

        # --- Save file ---
        savefile.save(model.state)
        prof.mark("Save")

        if kt % 24 == 0:
            print(prof.format_step(kt, cfg.nt))

    print(prof.format_summary())
    savefile.close()