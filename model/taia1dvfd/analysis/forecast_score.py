import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import butter, sosfiltfilt


def compute_forecast_score_uv(u_model_mean, u_model_var, z_model, t_model, n_ens,
                              u_obs, z_obs, n_obs, norm='first', filt=True, ddof=0,
                              return_terms=False):
    um = interp1d(z_model, u_model_mean)(z_obs)
    uv = interp1d(z_model, u_model_var)(z_obs)

    MSE = np.sum((um.real - u_obs.real)**2 + (um.imag - u_obs.imag)**2, axis=1) / (2 * n_obs)
    VAR = np.sum(uv.real + uv.imag, axis=1) / (2 * n_obs)
    return _EPS(MSE, VAR, n_ens, t_model, norm, filt, ddof, return_terms)


def compute_forecast_score_phi(phi_model_mean, phi_model_var, z_model, t_model, n_ens,
                               phi_obs, z_obs=None, interp=True,
                               norm='first', filt=True, ddof=0, return_terms=False):
    """phi_obs: (nt, n_obs) with interp, or (nt,) for a scalar such as the MLD."""
    if interp:
        phim = interp1d(z_model, phi_model_mean)(z_obs)
        phiv = interp1d(z_model, phi_model_var)(z_obs)
    else:
        phim, phiv = np.asarray(phi_model_mean), np.asarray(phi_model_var)

    phi_obs = np.asarray(phi_obs)
    if phim.shape != phi_obs.shape:
        raise ValueError(f'model {phim.shape} and obs {phi_obs.shape} shapes differ')

    MSE = (phim - phi_obs)**2
    VAR = phiv
    if MSE.ndim == 2:                       # several observation points
        MSE, VAR = MSE.mean(axis=1), VAR.mean(axis=1)
    return _EPS(MSE, VAR, n_ens, t_model, norm, filt, ddof, return_terms)


def _EPS(mse, var, n_ens, t_model, norm='first', filt=True, ddof=0, return_terms=False):
    
    cff = (n_ens + 1) / (n_ens - 1) if ddof == 0 else (n_ens + 1) / n_ens
    mse = np.asarray(mse, dtype=float)
    var_c = cff * np.asarray(var, dtype=float)
    eps = np.abs(mse - var_c)
    t_days = np.asarray(t_model) / 86400

    if filt:
        dt = t_days[1] - t_days[0]
        wn = 1.0 / (0.5 / dt)               # cutoff 1 day^-1 / Nyquist
        if not 0 < wn < 1:
            raise ValueError(f'output step {dt:.3g} d too coarse for a 1 d^-1 cutoff')
        sos = butter(5, wn, output='sos')
        eps = sosfiltfilt(sos, eps)
        if return_terms:
            mse, var_c = sosfiltfilt(sos, mse), sosfiltfilt(sos, var_c)

    if isinstance(norm, str) and norm == 'first':
        if eps[0] == 0:
            raise ValueError('eps(t0) = 0: use norm=None or an explicit normalisation')
        scale = eps[0]
    elif norm is not None and norm is not False:
        scale = norm
    else:
        scale = 1.0

    with np.errstate(divide='ignore', invalid='ignore'):
        eps = eps / scale

    if return_terms:
        return t_days, eps, {'mse': mse, 'var': var_c, 'scale': scale}
    return t_days, eps


def score_summary(t_days, eps, terms, t_start=10.0):
    """Time means over t >= t_start (days), for the score table."""
    idx = t_days >= t_start
    mse = terms['mse'][idx].mean()
    var = terms['var'][idx].mean()
    scale = terms['scale'] if np.isscalar(terms['scale']) else terms['scale'][idx].mean()
    return {
        'MSE':        mse,
        'VAR':        var,
        'VAR/MSE':    var / mse,              # < 1 under-dispersed, > 1 over-dispersed
        'signed':     (mse - var) / scale,    # > 0 under-, < 0 over-dispersed
        'eps_mean':   eps[idx].mean(),        # time mean of the plotted curve
    }