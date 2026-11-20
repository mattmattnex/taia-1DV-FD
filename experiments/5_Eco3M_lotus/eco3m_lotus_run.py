import torch
import numpy as np
import os
import netCDF4 as nc

from taia1dvfd.core  import Config, Grid, State, Model
from taia1dvfd.utils import SaveNetCDF, StepProfiler, Animator, Panel

import taia1dvfd.physics as phy
import taia1dvfd.mixing  as mix
import taia1dvfd.equation_of_state as eos
import taia1dvfd.biogeochem as bio
import taia1dvfd.forcing as forcing

ANIM = False
TEMP_FILE   = '/home/mnex/Desktop/1DV_FD/LOTUS/obs/lotusnstemp.epic'
SAVE_FOLDER = "/home/mnex/Desktop/1DV_FD/experiments/5_Eco3M_lotus/output/"
SEED_FORCING = 1234
SEED_NOISE   = 5678

# ======================================================================
def make_config() -> Config:
    cfg = Config()

    cfg.ne = 20
    cfg.nz = 64
    cfg.Hmin = 0.0
    cfg.Hmax = 150.0
    cfg.lat = 34.98
    cfg.device = "cuda:0"

    cfg.dt = 1800
    cfg.tend = 3600.0 * 24 * 160
    cfg.wind_Ta = 6 * 3600

    # Moderate, steady wind
    cfg.wind_mean = 4.0 + 0.0j
    cfg.wind_std = [0.5, 0.5]

    cfg.wave_mixing = True
    cfg.tke_wave_mix = False
    cfg.gls_wave_mix = False

    # On for the ensemble
    cfg.temp_stochastic = True
    cfg.salt_stochastic = True
    cfg.vel_stochastic  = True
    cfg.bgc_stochastic  = False

    cfg.Tsurf = 14.0
    cfg.Tbotm = 10.0
    cfg.Ssurf = 38.0
    cfg.Sbotm = 38.0

    cfg.tke_pdl = True
    cfg.tke_ebb = 67.83
    cfg.tke_azm_min = 1.2e-4
    cfg.tke_azt_min = 1.2e-5

    cfg.nbgc = 15

    cfg.save_folder = SAVE_FOLDER
    cfg.file_name = "eco3M_lotus.nc"
    cfg.var_to_save = [
        'mld', 'N2', 'temp',
        'bgc.phy_PHY_C', 'bgc.phy_PHY_Chl', 'bgc.phy_PHY_cell',
        'bgc.mid_DIN_N', 'bgc.mod_MODLN_N', 'bgc.mop_MOP_N',
        'bgc.bac_BAC_C', 'bgc.bac_BAC_N', 'bgc.zoo_COP_C', 'bgc.zoo_COP_N',
        'bgc.phy_PHY_N',
    ]
    cfg.finalize()
    return cfg

# ======================================================================
# Forcing: lotus
# ======================================================================
def initial_temperature(model: Model, z_profile: np.ndarray, T_profile: np.ndarray) -> torch.Tensor:
    """Observed profile (z <= 0, any order) -> cell centres, constant beyond the data."""
    zc = model.grid.z_c.reshape(-1).cpu().numpy()
    order = np.argsort(z_profile)
    T_c = np.interp(zc, z_profile[order], T_profile[order])
    T0 = torch.tensor(T_c, **model.config.farr_kwargs)
    return T0.reshape(1, -1, 1).expand(model.config.ne, -1, 1).clone()

cfg_ref   = make_config()
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

# ======================================================================
# Model init
# ======================================================================
cfg = make_config()
grid = Grid(cfg)
state = State(cfg)
model = Model(cfg, grid, state)
savefile = SaveNetCDF(model)
prof = StepProfiler(sync=True)

phy.noise_LU_init(model)
phy.wind_ornstein_init(model)
mix.tke_init(model)
mix.mld_init(model)

eos.temperature_equation_init(model, temp_init=initial_temperature(model, z_profile, T_profile))
eos.salinity_equation_init(model)
eos.eos_init(model)

Eco3M = bio.bgc_coupling_init(model)

savefile.init_file(model.state)
savefile.save(model.state)

# ======================================================================
# Main Loop 
# ======================================================================
for kt in range(cfg.nt):
    prof.start()

    # --- Wind ---
    phy.wind_ornstein_step(model)
    u_star = phy.compute_ustar(model).squeeze(-1)
    prof.mark("Wind")

    # --- Vertical Mixing ---
    mix.tke_step(u_star, model)
    prof.mark("Mix")
    
    # --- Tracers and EOS ---
    eos.temperature_equation_step(model, heat_flux_ns=Q_ns[:, kt], heat_flux_sw=Q_sw[kt])
    eos.salinity_equation_step(model)
    eos.eos_step(model)
    mix.mld_step(model)
    prof.mark("Eos")

    # --- Momentum ---
    phy.noise_LU_step(model)
    phy.momentum_equation_step(model)
    prof.mark("Mom")
    
    # --- Eco3M ---
    bio.bgc_coupling_step(model, Eco3M, Q_sw[kt].cpu().numpy())
    model.state.bgc[:] = eos.tracer_batch_equation_step(model, model.state.bgc, stochastic=cfg.bgc_stochastic)
    prof.mark("Eco3M")
    
    # --- End step ---
    model.state.time += model.config.dt
    savefile.save(model.state)
    
    if kt % 5 == 0:
        print(prof.format_step(kt, cfg.nt))
        
print(prof.format_summary())
savefile.close()