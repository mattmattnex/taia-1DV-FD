import torch

from taia1dvfd.core  import Config, Grid, State, Model, Operators
from taia1dvfd.utils import SaveNetCDF, StepProfiler

import taia1dvfd.physics as phy
import taia1dvfd.mixing  as mix
import taia1dvfd.equation_of_state as eos

for mod_name in ['k-kl', 'k-epsilon', 'k-omega', 'gen']:
    for nz_val in [512, 256]:
        for dt_val in [36, 360, 3600]:

            cfg = Config()

            cfg.ne   = 50
            cfg.nz   = nz_val
            cfg.Hmin = 1.0
            cfg.Hmax = 50.0
            cfg.lat  = 34.0
            cfg.device = "cuda:0"

            # --- Time ---
            cfg.dt   = dt_val
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

            # --- GLS model ---
            cfg.gls_model = mod_name
            cfg.gls_stabf = 'CA'
            cfg.gls_bc    = 'dirichlet'
            cfg.gls_azm_min = 1.2e-4
            cfg.gls_azt_min = 1.2e-5

            # --- Finalization ---
            cfg.save_folder = "/home/mnex/Desktop/1DV_FD/experiments/1_Kato_Phillips/output/"
            cfg.file_name = f"{cfg.gls_model}_dt{cfg.dt}_nz{cfg.nz}_ne{cfg.ne}.nc"
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

            mix.gls_init(model)
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
                model.state.ua[:] = 0.01 / cfg.sqrt_cdrho
                u_star = phy.compute_ustar(model).squeeze(-1)
                prof.mark("Wind")
                
                # --- Vertical Mixing ---
                mix.gls_step(u_star, model)
                prof.mark("GLS")
                
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