import math
import numpy as np
import torch
import netCDF4 as nc

LOTUS_MET_FILE = '/home/mnex/Desktop/1DV_FD/LOTUS/obs/lotusmet.epic'
LOTUS_BAD_VALUE = 1e30             # EPIC fill value is 9.96921e+36

def read_lotus_met(n_days: float = 160, path: str = LOTUS_MET_FILE) -> dict:
    """
    Read the LOTUS VAWR record (float64 numpy arrays at observation times).
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
