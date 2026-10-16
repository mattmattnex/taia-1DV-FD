import numpy as np
import scipy as scp

def compute_forecast_score_uv(u_model_mean, u_model_var, z_model, t_model, n_ens,
                              u_obs, z_obs, n_obs):
    
        # --- Interpolation at z_lotus ---
        um_int = scp.interpolate.interp1d(z_model, u_model_mean)(z_obs)
        uv_int = scp.interpolate.interp1d(z_model, u_model_var)(z_obs)
    
        # --- MSE and VAR ---
        MSE = (um_int.real - u_obs.real)**2 + (um_int.imag - u_obs.imag)**2
        VAR = uv_int.real + uv_int.imag

        MSE = np.sum(MSE, axis = 1) / (2 * n_obs)
        VAR = np.sum(VAR, axis = 1) / (2 * n_obs)

        # --- Forecast score ---
        cff = (n_ens + 1) / (n_ens-1)
        EPS = np.abs(MSE - cff * VAR)
        EPS /= EPS[0]

        # --- Filtering ---
        t_days = t_model/86400
        dt = t_days[1] - t_days[0] # days
        fs = 1./dt                 # samp. freq.
        fc = 1.                    # cutoff freq. (1/days)
        wid = fc/(fs/2)
        num_filt, denom_filt = scp.signal.butter(5, wid, 'low')

        return t_days, scp.signal.filtfilt(num_filt, denom_filt, EPS)