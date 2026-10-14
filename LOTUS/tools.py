import numpy as np  
import matplotlib.pyplot as plt
import netCDF4 as nc
import scipy as scp

class LotusTools:
    def __init__(self, output_folder):
        self.output_folder = output_folder
        
        # ------------ Observations ------------
        # Velocity from Price et al. 1987
        self.CWC_mean = np.array([0.046,  0.028,  0.020,  0.004]) #, np.nan,  0.006,  0.011])    
        self.CWC_std  = np.array([0.012,  0.007,  0.007,  0.004]) #, np.nan,  0.003,  0.006])
        self.DWC_mean = np.array([0.010, -0.003, -0.002, -0.005]) #, np.nan, -0.002,  0.007])
        self.DWC_std  = np.array([0.007,  0.004,  0.005,  0.004]) #, np.nan,  0.004,  0.006])
        z_lotus  = np.array([-5.0,   -10.0,  -15.0,  -25.0]) #,  -50.0,  -75.0, -100.0])
        # -------------------------------------

        self.uobs =  self.DWC_mean[np.newaxis, :]
        self.vobs = -self.CWC_mean[np.newaxis, :]
        self.nobs =  self.CWC_mean.shape[0]
        self.zobs = z_lotus
        
    def get_var(self, 
                file_name:str, 
                var_name: str, 
                cplx: bool = False, 
                kt:int | None = None):
        
        with nc.Dataset(file_name) as data:
            # --- Variables ---
            if cplx:
                var_m = data[f'{var_name}_mean_re'][:] + 1j*data[f'{var_name}_mean_im'][:]
                var_v = data[f'{var_name}_var_re'][:]  + 1j*data[f'{var_name}_var_im'][:]
            else:
                var_m = data[f'{var_name}_mean'][:]
                var_v = data[f'{var_name}_var'][:]

            # --- Coordinates ---
            z_c = data['zc'][:]
            z_w = data['zw'][:]
            t   = data['time'][:]

            # --- Find shape ---
            z = z_c if var_m.shape[1] == z_c.shape[0] else z_w
            
        if kt is not None:
            return z, var_m[kt, :], var_v[kt, :]
        else:
            return z, t, var_m, var_v
        
    def forecast_score(self, 
                       file_name: str, 
                       ne: int = 100):
        
        # --- Get var data ---
        z, t, um, uv = self.get_var(file_name, 'ue', cplx=True)
        
        # --- Interpolation at z_lotus ---
        um_int = scp.interpolate.interp1d(z, um)(self.zobs)
        uv_int = scp.interpolate.interp1d(z, uv)(self.zobs)
        
        # --- MSE and VAR ---
        MSE = (um_int.real - self.uobs)**2 + (um_int.imag - self.vobs)**2
        VAR = uv_int.real + uv_int.imag

        MSE = np.sum(MSE, axis = 1) / (2 * self.nobs)
        VAR = np.sum(VAR, axis = 1) / (2 * self.nobs)

        # --- Forecast score ---
        cff = (ne + 1) / (ne-1)
        EPS = np.abs(MSE - cff * VAR)
        EPS /= EPS[0]

        # --- Filtering ---
        t_days = t/86400
        dt = t_days[1] - t_days[0] # days
        fs = 1./dt                 # samp. freq.
        fc = 1.                    # cutoff freq. (1/days)
        wid = fc/(fs/2)
        num_filt, denom_filt = scp.signal.butter(5, wid, 'low')

        return t_days, scp.signal.filtfilt(num_filt, denom_filt, EPS)

    def plot_profiles(self, 
                      files_dict: dict, 
                      hmax: float = 50.0, 
                      figsize: tuple = (10, 4), 
                      colors = None):
        
        fig, ax = plt.subplots(1, 2, figsize = figsize, sharey = True)

        kplot = 0
        for label, file in files_dict.items():
            path = self.output_folder + file
            z, uem, uev = self.get_var(path, 'ue', kt = -1, cplx = True)
            
            if colors is not None:
                color = colors[kplot]
            else:
                color = f'C{kplot}'
            
            # u-profile
            ax[0].plot(uem.real, z, label = label, color = color)
            ax[0].fill_betweenx(z, uem.real - uev.real**0.5,
                                uem.real + uev.real**0.5,
                                alpha = 0.1, color = color)
            
            # v-profile
            ax[1].plot(uem.imag, z, label = label, color = color)
            ax[1].fill_betweenx(z, uem.imag - uev.imag**0.5,
                                uem.imag + uev.imag**0.5,
                                alpha = 0.1, color = color)
            
            kplot += 1
            
        # --- LOTUS obs errorbar ---
        ax[0].errorbar(self.uobs[0], self.zobs, xerr = self.DWC_std, color = 'k', ecolor = 'k', fmt = 'o', label = 'LOTUS', capsize=5)
        ax[0].set_title(r'$u_e(z)$')
        ax[0].set_xlabel(r'$[m.s^{-1}]$')
        ax[0].set_ylabel('z (m)')
        ax[0].set_xlim(-0.02, 0.04)
        ax[0].set_ylim(-hmax, 1)
        ax[0].grid()
        
        ax[1].errorbar(self.vobs[0], self.zobs, xerr = self.CWC_std, color = 'k', ecolor = 'k', fmt = 'o', label = 'LOTUS', capsize=5)
        ax[1].set_title(r'$v_e(z)$')
        ax[1].set_xlabel(r'$[m.s^{-1}]$')
        ax[1].set_xlim(-0.06, 0.01)
        ax[1].set_ylim(-hmax, 1)
        ax[1].grid()
        
        handles, labels = ax[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='center left', bbox_to_anchor=(1.0, 0.5), ncol=1)
        fig.tight_layout()
    
    def plot_forecast(self, 
                       files_dict: dict,
                       ne: int,
                       ylim: tuple = (-0.05, 0.3),
                       figsize: tuple = (10, 3),
                       colors = None):
        
        fig, ax = plt.subplots(1, 1, figsize = figsize)

        kplot = 0
        for label, file in files_dict.items():
            
            if colors is not None:
                color = colors[kplot]
            else:
                color = f'C{kplot}'
            
            path = self.output_folder + file
            t, eps = self.forecast_score(path, ne = ne)
            ax.plot(t, eps, label = label, color = color)
            kplot += 1
        
        ax.set_xlabel('Time (days)')
        ax.set_ylabel(r'$\epsilon / \epsilon_0$')
        ax.set_ylim(ylim)
        ax.grid()
        ax.set_title(r'Normalized forecast sore $\varepsilon / \varepsilon_0$')
        
        handles, labels = ax.get_legend_handles_labels()
        fig.legend(handles, labels, loc='lower center', ncol=len(labels), bbox_to_anchor=(0.5, -0.1))
        fig.tight_layout()
        
    def plot_mixing(self, 
                      files_dict: dict, 
                      hmax: float = 50.0, 
                      figsize: tuple = (4, 5),
                      colors = None):
        
        fig, ax = plt.subplots(1, 1, figsize = figsize, sharey = True)

        kplot = 0
        for label, file in files_dict.items():
            path = self.output_folder + file
            z, azm_m, azm_v = self.get_var(path, 'az_m', kt = -1, cplx = False)
            
            if colors is not None:
                color = colors[kplot]
            else:
                color = f'C{kplot}'
            
            # u-profile
            ax.plot(azm_m, z, label = label, color = color)
            ax.fill_betweenx(z, azm_m - azm_v**0.5,
                                azm_m + azm_v**0.5,
                                alpha = 0.1, color = color)
            
            kplot += 1
            
        ax.set_title(r'$a_{zz}^m(z)$')
        ax.set_xlabel(r'$[m^2.s^{-1}]$')
        ax.set_ylabel('z (m)')
        ax.set_ylim(-hmax, 1)
        ax.grid()

        handles, labels = ax.get_legend_handles_labels()
        fig.legend(handles, labels, loc='center left', bbox_to_anchor=(1.0, 0.5), ncol=1)
        fig.tight_layout()