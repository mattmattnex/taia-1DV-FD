
import os
import json
import numpy as np
import torch
import netCDF4 as nc

from taia1dvfd.core import Config, Grid, State, Model
import taia1dvfd.physics as phy
import taia1dvfd.mixing as mix
import taia1dvfd.equation_of_state as eos
import taia1dvfd.forcing as forcing

OUT       = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output')
TEMP_FILE = "/home/mnex/Desktop/1DV_FD/LOTUS/obs/lotusnstemp.epic"
CLOSURES  = ['k-kl', 'k-epsilon', 'tke']        # add 'k-omega', 'gen' if wanted
N_DAYS    = 40
NE        = 100
DTS       = [1800, 600, 300]
SEED_FORCING, SEED_NOISE = 1234, 5678
REC_HOURS = 3.0
T_MAX     = 50.0                                 # degC, above which a member has blown up

RUNS = {f'{s}_dt{dt}': dict(scheme=s, dt=dt, temp_noise=True) for s in ('lagged', 'current') for dt in DTS}
RUNS['noT_dt1800'] = dict(scheme='lagged', dt=1800, temp_noise=False)


def make_config(mix_model, dt, temp_noise):
    cfg = Config()
    cfg.ne, cfg.nz = NE, 300
    cfg.Hmin, cfg.Hmax = 0.0, 150.0
    cfg.lat = 34.98
    cfg.device = "cuda:0"
    cfg.dt = dt
    cfg.tend = 86400.0 * N_DAYS
    cfg.wind_Ta = 6 * 3600
    cfg.wave_mixing = cfg.tke_wave_mix = cfg.gls_wave_mix = True
    cfg.vel_stochastic = True
    cfg.temp_stochastic = temp_noise
    cfg.salt_stochastic = False
    cfg.temp_bc_top = cfg.temp_bc_btm = "neumann"
    cfg.Ssurf = cfg.Sbotm = 36.5
    cfg.gls_stabf = 'CA'
    cfg.gls_model = mix_model
    cfg.save_folder, cfg.file_name, cfg.var_to_save = OUT, 'unused.nc', []
    cfg.finalize()
    return cfg


def mixing(model, mix_model, u_star=None):
    kind = 'tke' if mix_model == 'tke' else 'gls'
    if u_star is None:
        getattr(mix, f'{kind}_init')(model)
    else:
        getattr(mix, f'{kind}_step')(u_star, model)


os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, 'runs.json'), 'w') as f:
    json.dump(dict(closures=CLOSURES, dts=DTS, runs=RUNS, n_days=N_DAYS, ne=NE,
                   rec_hours=REC_HOURS, t_max=T_MAX), f, indent=2)

with nc.Dataset(TEMP_FILE) as f:
    z_prof = -np.asarray(f.variables['depth'][:])
    T_prof = np.asarray(f.variables['temp'][0, :, 0, 0])

forcing_cache = {}

for mix_model in CLOSURES:
    for name, r in RUNS.items():
        print(f'===== {mix_model} | {name} =====')
        cfg = make_config(mix_model, r['dt'], r['temp_noise'])
        model = Model(cfg, Grid(cfg), State(cfg))
        st = model.state

        # --- Forcing: one per dt, generated with the same seed ---
        if r['dt'] not in forcing_cache:
            torch.manual_seed(SEED_FORCING)
            met = forcing.read_lotus_met(n_days=N_DAYS + 1)
            ua = forcing.generate_lotus_wind(model, met, window_len=24)
            forcing_cache[r['dt']] = (ua, *forcing.generate_lotus_heat_flux(model, met, ua))
        ua_f, Q_ns, Q_sw = forcing_cache[r['dt']]

        # --- Initial state ---
        torch.manual_seed(SEED_NOISE)
        phy.noise_LU_init(model)
        mixing(model, mix_model)
        mix.mld_init(model)
        zc = model.grid.z_c.reshape(-1).cpu().numpy()
        order = np.argsort(z_prof)
        T0 = torch.tensor(np.interp(zc, z_prof[order], T_prof[order]), **cfg.farr_kwargs)
        eos.temperature_equation_init(model, temp_init=T0.reshape(1, -1, 1).expand(NE, -1, 1).clone())
        eos.salinity_equation_init(model)
        eos.eos_init(model)

        def tracers(kt):
            eos.temperature_equation_step(model, heat_flux_ns=Q_ns[:, kt], heat_flux_sw=Q_sw[kt])
            eos.salinity_equation_step(model)
            eos.eos_step(model)

        rec_every = max(1, int(round(REC_HOURS * 3600 / cfg.dt)))
        rec = {k: [] for k in ('t', 'blow', 'mld_m', 'mld_s', 'sst')}

        # --- Main loop ---
        for kt in range(cfg.nt):
            st.ua[:, 0, 0] = ua_f[:, kt]
            u_star = phy.compute_ustar(model).squeeze(-1)
            mixing(model, mix_model, u_star)

            if r['scheme'] == 'lagged':
                tracers(kt)
                phy.noise_LU_step(model)
            else:
                phy.noise_LU_step(model)
                tracers(kt)

            mix.mld_step(model)
            phy.momentum_equation_step(model)
            st.time += cfg.dt

            if kt % rec_every == 0:
                T = st.temp[..., 0]
                bad = (~torch.isfinite(T) | (T.abs() > T_MAX)).any(dim=1)
                ok = ~bad
                mld = st.mld.reshape(-1).abs()[ok]
                rec['t'].append(float(st.time) / 86400)
                rec['blow'].append(bad.float().mean().item())
                rec['mld_m'].append(mld.mean().item() if ok.any() else np.nan)
                rec['mld_s'].append(mld.std().item() if ok.sum() > 1 else np.nan)
                rec['sst'].append(T[ok, 0].median().item() if ok.any() else np.nan)

        np.savez(os.path.join(OUT, f'{mix_model}_{name}.npz'), **{k: np.array(v) for k, v in rec.items()})
        print(f"   max blow {max(rec['blow']):.2f} | final MLD {rec['mld_m'][-1]:.1f} m")
        del model, st
        torch.cuda.empty_cache()
