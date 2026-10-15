import numpy as np
import torch
from taia1dvfd.core import Model

from .utils import _ornstein_anomaly

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