import torch

from taia1dvfd.core  import Config, Grid, State, Model
from taia1dvfd.utils import SaveNetCDF, StepProfiler

import taia1dvfd.physics as phy
import taia1dvfd.mixing  as mix
import taia1dvfd.equation_of_state as eos

# =====================================================================
#  Kato-Phillips experiment
#  Constant stress u* = 1 cm/s, no rotation, no buoyancy flux, no waves,
#  linear stratification N2 = 1e-4 s-2.  Reference: D(t) = 1.05 u* sqrt(t/N).
# =====================================================================

SEED      = 1234
MODES     = ['deterministic', 'stochastic']
CLOSURES  = ['k-kl', 'k-epsilon', 'k-omega', 'gen', 'tke']
NZ_VALUES = [512, 128, 32]
DT_VALUES = [36, 360, 3600]

HMAX      = 80.0       # deep enough that the bottom (artificial) stays far below the final MLD (~35 m)
OUTPUT_DT = 900.0      # output interval [s] (or dt if dt is larger) -> same output times for all dt
N2_INIT   = 1e-4       # initial stratification [s-2]
U_STAR    = 0.01       # friction velocity [m/s]

for mode in MODES:
    for closure in CLOSURES:
        for nz in NZ_VALUES:
            for dt in DT_VALUES:
                print(f'Running : {closure} with nz={nz}, dt={dt} s in mode {mode}')

                # --- Same random stream for every configuration ---
                torch.manual_seed(SEED)

                cfg = Config()
                cfg.nz   = nz
                cfg.Hmin = 0.0
                cfg.Hmax = HMAX
                cfg.lat  = 0.0                      # f = 0: no rotation
                cfg.device = "cuda:0"

                # --- Time ---
                cfg.dt   = dt
                cfg.tend = 3600.0 * 30
                cfg.save_every = max(1, int(round(OUTPUT_DT / dt)))

                # --- Stochastic bool ---
                cfg.temp_stochastic = False # No stochastic on tracers for T and S
                cfg.salt_stochastic = False
                                    
                if mode == 'deterministic':
                    cfg.ne = 1
                    cfg.vel_stochastic  = False
                elif mode == 'stochastic':
                    cfg.ne = 500
                    cfg.vel_stochastic  = True

                # --- No waves ---
                cfg.us0          = 0.0
                cfg.wave_mixing  = False
                cfg.tke_wave_mix = False
                cfg.gls_wave_mix = False

                # --- EOS: linear T profile giving N2 = N2_INIT, uniform S ---
                cfg.Tsurf = 13.0
                cfg.Tbotm = cfg.Tsurf - (N2_INIT / (cfg.g * cfg.beta_T)) * (cfg.Hmax - cfg.Hmin)
                cfg.Ssurf = 38.0
                cfg.Sbotm = 38.0

                # --- No surface buoyancy flux: zero-flux tracers ---
                cfg.temp_bc_top = cfg.temp_bc_btm = "neumann"
                cfg.salt_bc_top = cfg.salt_bc_btm = "neumann"

                # --- Mixing model ---
                if closure in mix._GLS_CLOSURES:
                    cfg.gls_model = closure
                    cfg.gls_stabf = 'CA'
                    cfg.gls_bc    = 'dirichlet'

                # --- Finalization ---
                cfg.save_folder = f"/home/mnex/Desktop/1DV_FD/experiments/1_Kato_Phillips/output/{mode}/"
                cfg.file_name   = f"{closure}_dt{cfg.dt}_nz{cfg.nz}_ne{cfg.ne}.nc"
                cfg.var_to_save = ['temp', 'mld', 'N2', 'tke', 'az_m', 'az_t']
                cfg.finalize()

                # --- Model objects ---
                grid      = Grid(cfg)
                state     = State(cfg)
                model     = Model(cfg, grid, state)
                savefile  = SaveNetCDF(model)
                prof      = StepProfiler()

                # --- Initializations ---
                phy.noise_LU_init(model)

                if closure in mix._GLS_CLOSURES:
                    mix.gls_init(model)
                elif closure == 'tke':
                    mix.tke_init(model)
                else:
                    raise ValueError(f"{closure} is not a valid closure")

                mix.mld_init(model)

                eos.temperature_equation_init(model)
                eos.salinity_equation_init(model)
                eos.eos_init(model)                     # computes rho and N2 from the initial T/S

                # --- Constant wind stress: u* = U_STAR ---
                model.state.ua[:] = U_STAR / cfg.sqrt_cdrho
                u_star = phy.compute_ustar(model).squeeze(-1)

                # --- Output file ---
                savefile.init_file(model.state)
                savefile.save(model.state)

                # ================== Main Loop ==================
                # Order : Mixing -> Noise -> T/S -> Momentum -> EOS -> MLD
                # ===============================================
                for kt in range(cfg.nt):
                    prof.start()

                    # --- Vertical Mixing ---
                    if closure in mix._GLS_CLOSURES:
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

                    # --- MLD ---
                    mix.mld_step(model)
                    prof.mark("mld")

                    # --- Update time ---
                    model.state.time += model.config.dt

                    # --- Save file ---
                    savefile.save(model.state)
                    prof.mark("Save")

                    if kt % 100 == 0:
                        print(prof.format_step(kt, cfg.nt))

                print(prof.format_summary())
                savefile.close()