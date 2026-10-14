import torch
torch.manual_seed(0)

from taia1dvfd.core  import Config, Grid, State, Model, Operators
from taia1dvfd.utils import SaveNetCDF, StepProfiler

import taia1dvfd.physics as phy
import taia1dvfd.mixing  as mix
import taia1dvfd.equation_of_state as eos

CLOSURES = ['k-kl', 'k-epsilon', 'k-omega', 'gen', 'tke']
NZ_VALUES = [500, 50, 25, 10]
DT_VALUES = [36, 360, 1200, 3600]
NENS = 500
# OUTPUT_FOLDER = 'deterministic'
OUTPUT_FOLDER = 'stochastic'

for closure in CLOSURES:
    for nz in NZ_VALUES:
        for dt in DT_VALUES:
            cfg = Config()

            cfg.ne   = NENS
            cfg.nz   = nz
            cfg.Hmin = 0.0
            cfg.Hmax = 50.0
            cfg.lat  = 0.0
            cfg.device = "cuda:0"

            # --- Time ---
            cfg.dt   = dt
            cfg.tend = 3600.0 * 30

            # --- Stochastic bool ---
            cfg.temp_stochastic = False
            cfg.salt_stochastic = False
            cfg.vel_stochastic  = True
            
            # --- No waves ---
            cfg.us0 = 0.0

            # --- EOS ---
            cfg.Tsurf  = 13
            cfg.Tbotm  = cfg.Tsurf - (1e-4 / (cfg.g * cfg.beta_T)) * (cfg.Hmax - cfg.Hmin)
            cfg.Ssurf = 38
            cfg.Sbotm = 38

            # --- Mixing model ---
            if closure in mix._GLS_CLOSURES:
                cfg.gls_model = closure
                cfg.gls_stabf = 'CA'
                cfg.gls_bc    = 'dirichlet'
                cfg.gls_azm_min = 1.2e-4
                cfg.gls_azt_min = 1.2e-5
            else:
                cfg.tke_azm_min = 1.2e-4
                cfg.tke_azt_min = 1.2e-5
                cfg.tke_wave_mix = False

            # --- Finalization ---
            cfg.save_folder = f"/home/mnex/Desktop/1DV_FD/experiments/1_Kato_Phillips/output/{OUTPUT_FOLDER}/"
            cfg.file_name = f"{closure}_dt{cfg.dt}_nz{cfg.nz}_ne{cfg.ne}.nc"
            cfg.var_to_save = ['temp', 'mld', 'N2']
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

            if closure in mix._GLS_CLOSURES:
                mix.gls_init(model)
            elif closure == 'tke':
                mix.tke_init(model)
            else:
                raise ValueError(f"{closure} is not a valid closure")

            mix.mld_init(model)

            eos.temperature_equation_init(model)
            eos.salinity_equation_init(model)
            eos.eos_init(model)

            savefile.init_file(model.state)
            savefile.save(model.state)

            # ================== Main Loop ==================
            # Order : Wind -> Mixing + Ops -> Noise -> T/S -> Momentum -> EOS -> MLD
            # ===============================================
            for kt in range(cfg.nt):
                prof.start()
                
                # --- Wind forcing ---
                model.state.ua[:] = 0.01 / cfg.sqrt_cdrho
                u_star = phy.compute_ustar(model).squeeze(-1)
                
                # --- Vertical Mixing ---
                if closure in mix._GLS_CLOSURES:
                    mix.gls_step(u_star, model)
                elif closure == 'tke':
                    mix.tke_step(u_star, model)
                operators.update_momentum_diffusion(state.az_m)
                operators.update_tracers_diffusion(state.az_t)
                prof.mark("Mix")
                
                # --- Noise ---
                phy.noise_LU_step(model)
                prof.mark("Sig")

                # --- Temp Sal ---
                eos.temperature_equation_step(model)
                eos.salinity_equation_step(model)
                prof.mark("T-S")
                
                # --- Momentum --- 
                phy.momentum_equation_step(model)
                prof.mark("Phy")

                # --- EOS ---
                eos.eos_step(model)

                # --- MLD ---
                mix.mld_step(model)
                prof.mark("mld")
                
                # --- Update time ---
                model.state.time += model.config.dt
                
                # --- Save file ---
                savefile.save(model.state)

                if kt%5 == 0:
                    print(prof.format_step(kt, cfg.nt))
            
            print(prof.format_summary())
            savefile.close()