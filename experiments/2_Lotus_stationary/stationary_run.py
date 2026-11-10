import torch

from taia1dvfd.core  import Config, Grid, State, Model
from taia1dvfd.utils import SaveNetCDF, StepProfiler

import taia1dvfd.physics as phy
import taia1dvfd.mixing  as mix
import taia1dvfd.equation_of_state as eos

SEED = 1234
CLOSURES = ['kpp', 'k-kl', 'k-epsilon', 'k-omega', 'gen', 'tke']
NE   = 1000
NZ   = 512
HMIN = 1.0
HMAX = 100.0
LAT  = 34.98
DT   = 3600.0
TEND = 3600.0 * 24 * 30

for closure in CLOSURES:
    print(f'Running : {closure}')

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
    cfg.wind_std  = [0.0, 1.0]
    
    # --- Stochastic bool ---
    cfg.vel_stochastic  = True
    cfg.temp_stochastic = True
    cfg.salt_stochastic = True

    # --- Waves ---
    cfg.wave_mixing  = True
    cfg.tke_wave_mix = False
    cfg.gls_wave_mix = False

    # --- EOS: neutral column (uniform T and S, N2 = 0) ---
    cfg.Tsurf = cfg.Tbotm = 20.0
    cfg.Ssurf = cfg.Sbotm = 38.0
    cfg.temp_bc_top = cfg.temp_bc_btm = "neumann"
    cfg.salt_bc_top = cfg.salt_bc_btm = "neumann"

    # --- Mixing model ---
    if closure == 'kpp':
        cfg.kpp_regularization = False
    elif closure in mix._GLS_CLOSURES:
        cfg.gls_model = closure
        cfg.gls_stabf = 'CA'
        cfg.gls_bc    = 'dirichlet'

    # --- Finalization ---
    cfg.save_folder = f"/home/mnex/Desktop/1DV_FD/experiments/2_Lotus_stationary/output/"
    cfg.file_name   = f"{closure}_ne{cfg.ne}.nc"
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

    if closure == 'kpp':
        mix.kpp_init(model)
    elif closure in mix._GLS_CLOSURES:
        mix.gls_init(model)
    elif closure == 'tke':
        mix.tke_init(model)
    else:
        raise ValueError(f"{closure} is not a valid closure")

    eos.temperature_equation_init(model)
    eos.salinity_equation_init(model)
    eos.eos_init(model)                 # N2 = 0 exactly (uniform anomalies)

    # --- Output file ---
    savefile.init_file(model.state)
    savefile.save(model.state)

    # ================== Main Loop ==================
    # Order : Wind -> Mixing -> Noise -> T/S -> Momentum -> EOS
    # ===============================================
    for kt in range(cfg.nt):
        prof.start()

        # --- Wind forcing ---
        phy.wind_ornstein_step(model)
        u_star = phy.compute_ustar(model).squeeze(-1)
        prof.mark("Wind")

        # --- Vertical Mixing ---
        if closure == 'kpp':
            mix.kpp_step(u_star, model)
        elif closure in mix._GLS_CLOSURES:
            mix.gls_step(u_star, model)
        else:
            mix.tke_step(u_star, model)
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