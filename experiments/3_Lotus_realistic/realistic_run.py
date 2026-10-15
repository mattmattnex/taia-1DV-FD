import numpy as np
import netCDF4 as nc
import torch

from taia1dvfd.core  import Config, Grid, State, Model
from taia1dvfd.utils import SaveNetCDF, StepProfiler

import taia1dvfd.physics as phy
import taia1dvfd.mixing  as mix
import taia1dvfd.equation_of_state as eos
import taia1dvfd.forcing as forcing

CLOSURES = [
    ('kpp', 'kpp'),
    ('gls', 'k-kl'),
    ('gls', 'k-epsilon'),
    ('gls', 'k-omega'),
    ('gls', 'gen'),
    ('tke', 'tke')
]

SEED_FORCING = 1234
SEED_NOISE   = 5678

def make_config(mix_model: str) -> Config:
    cfg = Config()

    cfg.ne   = 1000
    cfg.nz   = 300
    cfg.Hmin = 1.0
    cfg.Hmax = 150.0
    cfg.lat  = 34.98
    cfg.device = "cuda:0"

    # --- Time ---
    cfg.dt   = 720
    cfg.tend = 3600.0 * 24 * 160
    cfg.wind_Ta = 6 * 3600
    
    # --- Wave ---
    cfg.wave_mixing  = True
    cfg.tke_wave_mix = True
    cfg.gls_wave_mix = True
    
    # --- Stochastic bool ---
    cfg.temp_stochastic = True
    cfg.salt_stochastic = False
    cfg.vel_stochastic  = True

    # --- Temperature BCs: heat flux at the surface, no flux at the bottom ---
    cfg.temp_bc_top = "neumann"
    cfg.temp_bc_btm = "neumann"

    # --- EOS (uniform salinity, Sargasso value) ---
    cfg.Ssurf = 36.5
    cfg.Sbotm = 36.5

    # --- GLS ---
    cfg.gls_stabf = 'CA'
    cfg.gls_model = mix_model

    # --- Output ---
    cfg.save_folder = "/home/mnex/Desktop/1DV_FD/experiments/3_Lotus_realistic/output/"
    cfg.file_name = f'lotus_stratif_{mix_model}_ne{cfg.ne}.nc'
    cfg.var_to_save = ['ue', 'temp', 'mld']
    cfg.finalize()
    return cfg

def initial_temperature(model: Model, z_profile: np.ndarray, T_profile: np.ndarray) -> torch.Tensor:
    """Observed profile (z <= 0, any order) -> cell centres, constant beyond the data."""
    zc = model.grid.z_c.reshape(-1).cpu().numpy()
    order = np.argsort(z_profile)
    T_c = np.interp(zc, z_profile[order], T_profile[order])
    T0 = torch.tensor(T_c, **model.config.farr_kwargs)
    return T0.reshape(1, -1, 1).expand(model.config.ne, -1, 1).clone()

# =====================================================================
#  Forcing identical for every closure
# =====================================================================
cfg_ref = make_config('kpp')
model_ref = Model(cfg_ref, Grid(cfg_ref), State(cfg_ref))

torch.manual_seed(SEED_FORCING)
met = forcing.read_lotus_met(n_days=cfg_ref.tend / 86400 + 1)
ua_forcing = forcing.generate_lotus_wind(model_ref, met, window_len=24)
Q_ns, Q_sw = forcing.generate_lotus_heat_flux(model_ref, met, ua_forcing)

# --- Initial temperature (13 May 1982) ---
with nc.Dataset('/home/mnex/Desktop/1DV_FD/LOTUS/obs/lotusnstemp.epic', 'r') as f:
    z_profile = -f.variables['depth'][:]
    T_profile = f.variables['temp'][0, :, 0, 0]

del model_ref

# =====================================================================
#  Closures
# =====================================================================
for prefix, mix_model in CLOSURES:
    cfg = make_config(mix_model)

    # --- Model objects ---
    grid   = Grid(cfg)
    state  = State(cfg)
    model    = Model(cfg, grid, state)
    savefile = SaveNetCDF(model)
    prof = StepProfiler()

    # --- Same LU noise sequence for every closure ---
    torch.manual_seed(SEED_NOISE)

    # --- Initializations ---
    phy.noise_LU_init(model)

    if mix_model == 'kpp':
        mix.kpp_init(model)
    elif mix_model == 'tke':
        mix.tke_init(model)
    elif mix_model in ['k-kl', 'k-epsilon', 'k-omega', 'gen']:
        mix.gls_init(model)
    else:
        raise ValueError('Mix model not correct')

    mix.mld_init(model)

    eos.temperature_equation_init(model, temp_init=initial_temperature(model, z_profile, T_profile))
    eos.salinity_equation_init(model)
    eos.eos_init(model)

    savefile.init_file(model.state)
    savefile.save(model.state)

    # ================== Main Loop ==================
    for kt in range(cfg.nt):
        prof.start()

        # --- Wind ---
        model.state.ua[:, 0, 0] = ua_forcing[:, kt]
        u_star = phy.compute_ustar(model).squeeze(-1)
        prof.mark("Wind")

        # --- Vertical Mixing ---
        if mix_model == 'kpp':
            mix.kpp_step(u_star, model)
        elif mix_model == 'tke':
            mix.tke_step(u_star, model)
        else:
            mix.gls_step(u_star, model)
        prof.mark("Mix")

        # --- Tracers and EOS ---
        eos.temperature_equation_step(model, heat_flux_ns=Q_ns[:, kt], heat_flux_sw=Q_sw[kt])
        eos.salinity_equation_step(model)
        eos.eos_step(model)
        prof.mark("eos")

        # --- MLD ---
        mix.mld_step(model)
        prof.mark("mld")

        # --- Physics and stochastic ---
        phy.noise_LU_step(model)
        phy.momentum_equation_step(model)
        prof.mark("Phy")

        # --- Update time ---
        model.state.time += model.config.dt

        # --- Save file ---
        savefile.save(model.state)

        if kt % 24 == 0:
            print(prof.format_step(kt, cfg.nt))

    print(prof.format_summary())
    savefile.close()
