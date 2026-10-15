import math
import numpy as np
import torch
import netCDF4 as nc
from taia1dvfd.core import Model

LOTUS_MET_FILE = '/home/mnex/Desktop/1DV_FD/LOTUS/obs/lotusmet.epic'
LOTUS_BAD_VALUE = 1e30             # EPIC fill value is 9.96921e+36

# --- Physical constants ---
STEFAN_BOLTZMANN = 5.67e-8         # W m-2 K-4
R_DRY_AIR = 287.05                 # J kg-1 K-1
CP_AIR = 1004.67                   # J kg-1 K-1
L_VAP = 2.5e6                      # J kg-1
KELVIN = 273.15


# =====================================================================
#  Observations
# =====================================================================
def read_lotus_met(n_days: float = 160, path: str = LOTUS_MET_FILE) -> dict:
    """
    Read the LOTUS VAWR record (float64 numpy arrays at observation times).

      t     : seconds since the first record, from the EPIC time / time2 pair
              (so gaps and irregular steps are handled, no fixed 900 s assumed)
      wnde, wndn (m/s), atmp, stmp (degC), srad (W/m2), bpr (hPa)

    Fill values are replaced by linear interpolation in time.
    """
    with nc.Dataset(path, 'r') as f:
        f.set_auto_mask(False)
        jday = f.variables['time'][:].astype(np.float64)       # true Julian day
        msec = f.variables['time2'][:].astype(np.float64)      # ms since 00:00 GMT
        t_all = (jday - jday[0]) * 86400 + (msec - msec[0]) * 1e-3

        n_obs = min(np.searchsorted(t_all, n_days * 86400) + 2, len(t_all))
        met = {'t': t_all[:n_obs]}
        for name in ['wnde', 'wndn', 'atmp', 'stmp', 'srad', 'bpr']:
            raw = f.variables[name][:n_obs, 0, 0, 0].astype(np.float64)
            met[name] = _fill_gaps(met['t'], raw, name)

    met['srad'] = np.maximum(met['srad'], 0.0)                 # no negative night-time offsets
    return met


def _fill_gaps(t: np.ndarray, x: np.ndarray, name: str) -> np.ndarray:
    bad = ~np.isfinite(x) | (np.abs(x) > LOTUS_BAD_VALUE)
    if bad.all():
        raise ValueError(f'{name}: no valid data')
    if bad.any():
        print(f'{name}: {bad.sum()} fill values interpolated')
        x = x.copy()
        x[bad] = np.interp(t[bad], t[~bad], x[~bad])
    return x


def _window_mean(t_obs: np.ndarray, x: np.ndarray, t_start: np.ndarray, dt: float) -> np.ndarray:
    """Mean of x over [t_start, t_start + dt] (trapezoidal cumulative integral)."""
    x_int = np.concatenate(([0.0], np.cumsum(0.5 * (x[1:] + x[:-1]) * np.diff(t_obs))))
    return (np.interp(t_start + dt, t_obs, x_int) - np.interp(t_start, t_obs, x_int)) / dt


def _ornstein_anomaly(ne: int, nt: int, dt: float, Ta: float, std: torch.Tensor,
                      fkw: dict) -> torch.Tensor:
    """Zero-mean OU anomaly (ne, nt), eta(0) = 0, stationary std = std (scalar or (nt,))."""
    c1 = math.exp(-dt / Ta)
    c2 = math.sqrt(1.0 - c1**2)
    std = torch.as_tensor(std, **fkw).expand(nt)

    xi = torch.randn(ne, nt, **fkw)
    eta = torch.zeros(ne, nt, **fkw)
    for kt in range(1, nt):
        eta[:, kt] = c1 * eta[:, kt-1] + c2 * std[kt] * xi[:, kt]
    return eta


# =====================================================================
#  Wind
# =====================================================================
def generate_lotus_wind(model: Model, met: dict, window_len: int = 24) -> torch.Tensor:
    """
    Random 10 m wind (ne, nt) complex:

      ua = m_obs(t) + eta_u + i eta_v

    m_obs  : block means of the observations (window_len records per block)
    eta_*  : zero-mean OU anomalies, memory cfg.wind_Ta, std = block std of each component

    Changes from the previous version:
      - the OU acts on the anomaly, so the ensemble mean follows the observations
        (relaxing x towards m(t) made the mean a lagged, low-passed copy of m)
      - one real N(0, 1) per component (complex randn has variance 1/2 per part, and
        (std_u + i std_v) * complex noise mixed the two components)
    """
    cfg = model.config
    fkw = cfg.farr_kwargs
    ckw = cfg.carr_kwargs
    ne, nt, dt = cfg.ne, cfg.nt, cfg.dt

    # --- Block statistics ---
    t_obs = met['t']
    n_blocks = len(t_obs) // window_len
    t_block = t_obs[:n_blocks * window_len].reshape(n_blocks, window_len)[:, 0]

    def block_stats(x):
        blocks = x[:n_blocks * window_len].reshape(n_blocks, window_len)
        return blocks.mean(axis=1), blocks.std(axis=1, ddof=1)

    ua_mean, ua_std = block_stats(met['wnde'])
    va_mean, va_std = block_stats(met['wndn'])

    # --- Block index of each model step ---
    t_model = np.arange(nt) * dt
    kb = np.clip(np.searchsorted(t_block, t_model, side='right') - 1, 0, n_blocks - 1)

    wind_mean = torch.tensor(ua_mean[kb] + 1j * va_mean[kb], **ckw)          # (nt,)
    eta_u = _ornstein_anomaly(ne, nt, dt, cfg.wind_Ta, torch.tensor(ua_std[kb], **fkw), fkw)
    eta_v = _ornstein_anomaly(ne, nt, dt, cfg.wind_Ta, torch.tensor(va_std[kb], **fkw), fkw)

    return wind_mean[None, :] + eta_u + 1j * eta_v


# =====================================================================
#  Heat flux (bulk formulae, option A: no external data)
# =====================================================================
def generate_lotus_heat_flux(model: Model, met: dict, ua_forcing: torch.Tensor,
                             rh_mean: float = 0.78,
                             rh_std: float = 0.05,
                             rh_Ta: float = 2 * 86400.,
                             cloud: float | None = None,
                             lw_lambda: float = 0.69,
                             ce: float = 1.15e-3,
                             ch_unstable: float = 1.13e-3,
                             ch_stable: float = 0.66e-3,
                             albedo: float = 0.06,
                             emissivity: float = 0.98,
                             wind_min: float = 0.5,
                             ):
    """
    Surface heat fluxes at the model steps, positive INTO the ocean (W/m2).

      Q_ns (ne, nt) : non-solar = -(Q_E + Q_H + Q_LW), per member (wind and humidity)
      Q_sw (nt,)    : net shortwave (1 - albedo) srad, deterministic, absorbed over depth

    Observed (window-averaged over each model step): T_a, T_s, srad, p.
    Wind speed: |ua_forcing| of each member.

    Missing in the VAWR record, replaced by:
      - humidity  : RH = rh_mean + OU anomaly (rh_std, memory rh_Ta), clipped to [0.3, 1]
                    -> the latent flux is the main forcing uncertainty (+-5 % RH ~ +-20 W/m2)
      - longwave  : Clark et al. (1974)
                    Q_LW = eps sig Ts^4 (0.39 - 0.05 sqrt(e_a)) (1 - lambda C^2) + 4 eps sig Ts^3 (Ts - Ta)
      - clouds C  : constant `cloud`, or from the daily shortwave relative to a clear-sky
                    envelope (31-day running max of daily means), C = (1 - Q/Q_clear) / 0.62

    Turbulent fluxes: constant transfer coefficients of Large & Pond (1982),
      Q_E = rho_a L_v C_E U (q_s - q_a),     q_s = 0.98 q_sat(T_s)
      Q_H = rho_a c_p C_H U (T_s - T_a),     C_H unstable (T_s > T_a) / stable
    Fluxes use the OBSERVED SST, so every closure receives exactly the same forcing.
    """
    cfg = model.config
    fkw = cfg.farr_kwargs
    ne, nt, dt = cfg.ne, cfg.nt, cfg.dt

    # --- Observations averaged over each model step [kt dt, (kt+1) dt] ---
    t_model = np.arange(nt) * dt
    t_obs = met['t']
    Ta   = torch.tensor(_window_mean(t_obs, met['atmp'], t_model, dt), **fkw)   # (nt,)
    Ts   = torch.tensor(_window_mean(t_obs, met['stmp'], t_model, dt), **fkw)
    srad = torch.tensor(_window_mean(t_obs, met['srad'], t_model, dt), **fkw)
    pres = torch.tensor(_window_mean(t_obs, met['bpr'],  t_model, dt), **fkw)

    # --- Cloud fraction ---
    if cloud is not None:
        cloud_frac = torch.full((nt,), cloud, **fkw)
    else:
        cloud_frac = torch.tensor(_cloud_from_shortwave(met, t_model), **fkw)

    # --- Member fields (ne, nt) ---
    wind_speed = ua_forcing.abs().to(fkw['dtype']).clamp(min=wind_min)
    rel_hum = (rh_mean + _ornstein_anomaly(ne, nt, dt, rh_Ta, rh_std, fkw)).clamp(0.3, 1.0)

    # --- Humidity ---
    e_sat_s = _sat_vapour_pressure(Ts)                    # hPa
    e_air   = rel_hum * _sat_vapour_pressure(Ta)
    q_surf  = 0.98 * _specific_humidity(e_sat_s, pres)    # 0.98: salinity effect
    q_air   = _specific_humidity(e_air, pres)
    rho_air = pres * 100 / (R_DRY_AIR * (Ta + KELVIN) * (1 + 0.61 * q_air))

    # --- Turbulent fluxes (positive upward = ocean loss) ---
    dT_sa = Ts - Ta
    ch = ch_stable + (ch_unstable - ch_stable) * (dT_sa > 0).to(dT_sa.dtype)
    Q_E = rho_air * L_VAP * ce * wind_speed * (q_surf - q_air)
    Q_H = rho_air * CP_AIR * ch * wind_speed * dT_sa

    # --- Net longwave (positive upward) ---
    Ts_K = Ts + KELVIN
    Q_LW = (emissivity * STEFAN_BOLTZMANN * Ts_K**4
            * (0.39 - 0.05 * torch.sqrt(e_air)) * (1 - lw_lambda * cloud_frac**2)
            + 4 * emissivity * STEFAN_BOLTZMANN * Ts_K**3 * dT_sa)

    # --- Fluxes into the ocean ---
    Q_ns = -(Q_E + Q_H + Q_LW)
    Q_sw = (1 - albedo) * srad

    print(f'Heat flux (mean over run): Q_sw = {Q_sw.mean():.0f}, Q_E = {-Q_E.mean():.0f}, '
          f'Q_H = {-Q_H.mean():.0f}, Q_LW = {-Q_LW.mean():.0f}, '
          f'net = {(Q_sw.mean() + Q_ns.mean()):.0f} W/m2')

    return Q_ns, Q_sw


def _sat_vapour_pressure(T: torch.Tensor) -> torch.Tensor:
    """Buck (1981), T in degC, returns hPa."""
    return 6.1121 * torch.exp(17.502 * T / (240.97 + T))


def _specific_humidity(e: torch.Tensor, p: torch.Tensor) -> torch.Tensor:
    return 0.622 * e / (p - 0.378 * e)


def _cloud_from_shortwave(met: dict, t_model: np.ndarray) -> np.ndarray:
    """Daily cloud fraction from measured shortwave vs a 31-day clear-sky envelope."""
    day_obs = (met['t'] // 86400).astype(int)
    n_days = day_obs.max() + 1
    sw_daily = np.bincount(day_obs, weights=met['srad'], minlength=n_days) \
             / np.maximum(np.bincount(day_obs, minlength=n_days), 1)

    half_window = 15
    sw_clear = np.array([sw_daily[max(0, d - half_window):d + half_window + 1].max()
                         for d in range(n_days)])
    cloud_daily = np.clip((1 - sw_daily / np.maximum(sw_clear, 1.0)) / 0.62, 0.0, 1.0)

    day_model = np.clip((t_model // 86400).astype(int), 0, n_days - 1)
    return cloud_daily[day_model]
