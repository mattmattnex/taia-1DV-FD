import numpy as np
import netCDF4 as nc
import torch
import os

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

RUN_TYPES = {
    'sto': dict(vel_stochastic=True,  temp_stochastic=True),
    'det': dict(vel_stochastic=False, temp_stochastic=False),
}

SEED_FORCING = 1234
SEED_NOISE   = 5678

SAVE_FOLDER = "/home/mnex/Desktop/1DV_FD/experiments/3_Lotus_realistic/output/"
TEMP_FILE   = '/home/mnex/Desktop/1DV_FD/LOTUS/obs/lotusnstemp.epic'


def make_config(mix_model: str, run_type: str = 'sto') -> Config:
    cfg = Config()

    cfg.ne   = 1000
    cfg.nz   = 300
    cfg.Hmin = 0.0
    cfg.Hmax = 150.0
    cfg.lat  = 34.98
    cfg.device = "cuda:0"

    # --- Time ---
    cfg.dt   = 1800
    cfg.tend = 3600.0 * 24 * 160
    cfg.wind_Ta = 6 * 3600

    # --- Wave ---
    cfg.wave_mixing    = True
    cfg.tke_wave_mix   = True
    cfg.gls_wave_mix   = True

    # --- Stochastic bool ---
    cfg.salt_stochastic = False
    for name, value in RUN_TYPES[run_type].items():
        setattr(cfg, name, value)

    # --- Temperature BCs: heat flux at the surface, no flux at the bottom ---
    cfg.temp_bc_top = "neumann"
    cfg.temp_bc_btm = "neumann"

    # --- EOS (Sargasso value) ---
    cfg.Ssurf = 36.5
    cfg.Sbotm = 36.5

    # --- GLS ---
    cfg.gls_stabf = 'CA'
    cfg.gls_model = mix_model

    # --- Output ---
    cfg.save_folder = SAVE_FOLDER
    cfg.file_name   = f'lotus_stratif_{mix_model}_{run_type}_ne{cfg.ne}.nc'
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


def mixing_init(model: Model, mix_model: str):
    if mix_model == 'kpp':
        mix.kpp_init(model)
    elif mix_model == 'tke':
        mix.tke_init(model)
    elif mix_model in ['k-kl', 'k-epsilon', 'k-omega', 'gen']:
        mix.gls_init(model)
    else:
        raise ValueError(f'Mix model not correct: {mix_model}')


def mixing_step(model: Model, mix_model: str, u_star: torch.Tensor):
    if mix_model == 'kpp':
        mix.kpp_step(u_star, model)
    elif mix_model == 'tke':
        mix.tke_step(u_star, model)
    else:
        mix.gls_step(u_star, model)


# =====================================================================
#  Forcing and initial state, identical for every closure and run type
# =====================================================================
cfg_ref   = make_config('kpp')
model_ref = Model(cfg_ref, Grid(cfg_ref), State(cfg_ref))

torch.manual_seed(SEED_FORCING)
met = forcing.read_lotus_met(n_days=cfg_ref.tend / 86400 + 1)
ua_forcing = forcing.generate_lotus_wind(model_ref, met, window_len=24)
Q_ns, Q_sw = forcing.generate_lotus_heat_flux(model_ref, met, ua_forcing)

# Saving forcing for analysis
np.savez(os.path.join(cfg_ref.save_folder, 'lotus_forcing.npz'),
    ua   = ua_forcing.cpu().numpy(),
    Q_ns = Q_ns.cpu().numpy(),
    Q_sw = Q_sw.cpu().numpy(),
    dt   = cfg_ref.dt,
    dz   = model_ref.grid.dz.cpu().numpy(),
    tend = cfg_ref.tend
)

# --- Initial temperature (13 May 1982) ---
with nc.Dataset(TEMP_FILE, 'r') as f:
    z_profile = -f.variables['depth'][:]
    T_profile = f.variables['temp'][0, :, 0, 0]

del model_ref


# =====================================================================
#  Runs: every closure, stochastic and deterministic
# =====================================================================
for run_type in RUN_TYPES:
    for prefix, mix_model in CLOSURES:
        print(f'\n===== {mix_model} | {run_type} =====')
        cfg = make_config(mix_model, run_type)

        # --- Model objects ---
        grid     = Grid(cfg)
        state    = State(cfg)
        model    = Model(cfg, grid, state)
        savefile = SaveNetCDF(model)
        prof     = StepProfiler()

        # --- Same random sequence for every closure and run type ---
        torch.manual_seed(SEED_NOISE)

        # --- Initializations ---
        phy.noise_LU_init(model)
        mixing_init(model, mix_model)
        mix.mld_init(model)

        eos.temperature_equation_init(model, temp_init=initial_temperature(model, z_profile, T_profile))
        eos.salinity_equation_init(model)
        eos.eos_init(model)

        savefile.init_file(model.state)
        savefile.save(model.state)

        # ================== Main Loop ==================
        for kt in range(cfg.nt):
            prof.start()

            # --- Wind and waves ---
            model.state.ua[:, 0, 0] = ua_forcing[:, kt]
            u_star = phy.compute_ustar(model).squeeze(-1)
            prof.mark("Wind")

            # --- Vertical Mixing ---
            mixing_step(model, mix_model, u_star)
            prof.mark("Mix")

            # --- Noise ---
            phy.noise_LU_step(model)
            prof.mark("Noise")

            # --- Tracers and EOS ---
            eos.temperature_equation_step(model, heat_flux_ns=Q_ns[:, kt], heat_flux_sw=Q_sw[kt])
            eos.salinity_equation_step(model)
            eos.eos_step(model)
            prof.mark("eos")

            # --- MLD ---
            mix.mld_step(model)
            prof.mark("mld")

            # --- Momentum ---
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

        # --- Free GPU memory before the next run ---
        del model, state, grid, savefile
        torch.cuda.empty_cache()