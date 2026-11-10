import numpy as np
import netCDF4 as nc

G       = 9.81
ALPHA_T = 2e-4      # must be cfg.beta_T of the model
N2_MIN  = 2e-5      # s-2, ~0.1 degC over 10 m with alpha = 2e-4
BAD     = 1e30      # EPIC fill value 9.96921e+36


# =====================================================================
#  Reading
# =====================================================================
def epic_days(f) -> np.ndarray:
    """Absolute EPIC time in days (true Julian day + ms since 00:00 GMT)."""
    return f.variables['time'][:].astype(np.float64) + f.variables['time2'][:].astype(np.float64) / 86400e3

def met_start_day(met_file: str) -> float:
    """Model t = 0 is the first record of the meteorological file (read_lotus_met)."""
    with nc.Dataset(met_file) as f:
        f.set_auto_mask(False)
        return epic_days(f)[0]
    
def _fill_missing_in_depth(depth, field):
    """Linear interpolation in depth of NaN sensors, profile by profile (no extrapolation)."""
    for kt in np.where(np.isnan(field).any(axis=1))[0]:
        valid = ~np.isnan(field[kt])
        if valid.sum() >= 2:
            inside = (depth >= depth[valid][0]) & (depth <= depth[valid][-1])
            field[kt, inside] = np.interp(depth[inside], depth[valid], field[kt, valid])
    return field


def read_lotus_temp(temp_file: str, t0_day: float, n_days: float):
    """
    Returns
      t_days : (nt,)       days since the model start
      depth  : (nsens,)    sensor depths, positive downward, increasing
      temp   : (nt, nsens) degC, missing sensors filled by interpolation in depth
    """
    with nc.Dataset(temp_file) as f:
        f.set_auto_mask(False)
        t_days = epic_days(f) - t0_day
        depth  = np.abs(f.variables['depth'][:].astype(np.float64))
        temp   = f.variables['temp'][:, :, 0, 0].astype(np.float64)

    order = np.argsort(depth)
    depth, temp = depth[order], temp[:, order]
    temp[~np.isfinite(temp) | (np.abs(temp) > BAD)] = np.nan

    keep = (t_days >= 0.0) & (t_days <= n_days)
    t_days, temp = t_days[keep], temp[keep]

    # --- Missing sensors: linear interpolation in depth (no extrapolation) ---
    temp = _fill_missing_in_depth(depth, temp)

    print(f'{np.isnan(temp).any(axis=1).sum()} profiles still incomplete after filling')
    return t_days, depth, temp

def read_lotus_vel(vel_file: str, t0_day: float, n_days: float):
    """
    Returns
      t_days : (nt,)       days since the model start
      depth  : (nsens,)    sensor depths, positive downward, increasing
      u      : (nt, nsens) eastward velocity  [m/s], missing sensors filled by interpolation in depth
      v      : (nt, nsens) northward velocity [m/s], idem
    """
    with nc.Dataset(vel_file) as f:
        f.set_auto_mask(False)
        t_days = epic_days(f) - t0_day
        depth  = np.abs(f.variables['depth'][:].astype(np.float64))
        u      = f.variables['east' ][:, :, 0, 0].astype(np.float64)
        v      = f.variables['north'][:, :, 0, 0].astype(np.float64)

    order = np.argsort(depth)
    depth, u, v = depth[order], u[:, order], v[:, order]

    # --- Bad values: a sensor is invalid if either component is bad ---
    bad = (~np.isfinite(u) | (np.abs(u) > BAD) |
           ~np.isfinite(v) | (np.abs(v) > BAD))
    u[bad] = np.nan
    v[bad] = np.nan

    # --- cm/s -> m/s ---
    u *= 1e-2
    v *= 1e-2

    keep = (t_days >= 0.0) & (t_days <= n_days)
    t_days, u, v = t_days[keep], u[keep], v[keep]

    u = _fill_missing_in_depth(depth, u)
    v = _fill_missing_in_depth(depth, v)

    print(f'{np.isnan(u).any(axis=1).sum()} velocity profiles still incomplete after filling')
    return t_days, depth, u, v

# =====================================================================
#  MLD on the sensor grid
# =====================================================================
def mld_sensor(temp: np.ndarray, depth: np.ndarray,
               alpha_T: float = ALPHA_T, N2_min: float = N2_MIN) -> np.ndarray:
    """
    temp  : (..., nsens) at depth (nsens,), positive downward, increasing
    returns MLD (...) in m, positive downward (NaN only if the profile has NaNs)
    """
    dz      = depth[1:] - depth[:-1]
    depth_m = 0.5 * (depth[1:] + depth[:-1])                     # interfaces between sensors

    N2 = G * alpha_T * (temp[..., :-1] - temp[..., 1:]) / dz     # > 0 when stable
    N2_valid = np.where(np.isnan(N2), -np.inf, N2)

    idx_max = np.argmax(N2_valid, axis=-1)
    N2_max  = np.take_along_axis(N2_valid, idx_max[..., None], axis=-1)[..., 0]

    mld = np.where(N2_max >= N2_min, depth_m[idx_max], depth[-1])
    return np.where(np.isnan(temp).any(axis=-1), np.nan, mld)


# =====================================================================
#  Common time averaging
# =====================================================================
def bin_average(t: np.ndarray, x: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """NaN-aware mean of x(t) in each bin [edges[i], edges[i+1])."""
    idx = np.digitize(t, edges) - 1
    out = np.full(len(edges) - 1, np.nan)
    for i in range(len(edges) - 1):
        sel = (idx == i) & np.isfinite(x)
        if sel.any():
            out[i] = x[sel].mean()
    return out