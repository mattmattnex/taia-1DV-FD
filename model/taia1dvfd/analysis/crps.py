import numpy as np
from scipy.stats import norm

def crps_gaussian(mu, var, y):
    """CRPS of the distribution N(mu, var) for the reference y.
    For var = 0 (deterministic run) it equals |y - mu|."""
    sig = np.sqrt(np.maximum(var, 0.0))
    err = y - mu
    with np.errstate(divide='ignore', invalid='ignore'):
        z = err / sig
        c = sig * (z * (2 * norm.cdf(z) - 1) + 2 * norm.pdf(z) - 1 / np.sqrt(np.pi))
    return np.where(sig > 0, c, np.abs(err))
 
def kp_metrics(t, m, v, D, t_min=2.0):
    """Relative metrics (fractions of D_m). Time means over t >= t_min."""
    sel  = t >= t_min
    bias = (m - D) / D
    crps = crps_gaussian(m, v, D) / D
    return dict(t=t, bias_t=bias, crps_t=crps, sprd_t=np.sqrt(v) / D,
                crps=crps[sel].mean(),              # time-mean relative CRPS
                bias_tf=bias[-1],                   # relative bias at t_f
                spread_tf=np.sqrt(v[-1]) / D[-1])   # relative spread at t_f
