# import numpy as np  
# import matplotlib.pyplot as plt
# import netCDF4 as nc
# import scipy as scp

# class LotusTools:
#     def __init__(self, output_folder):
#         self.output_folder = output_folder
        
#         # ------------ Observations ------------
#         # Velocity from Price et al. 1987
#         self.CWC_mean = np.array([0.046,  0.028,  0.020,  0.004]) #, np.nan,  0.006,  0.011])    
#         self.CWC_std  = np.array([0.012,  0.007,  0.007,  0.004]) #, np.nan,  0.003,  0.006])
#         self.DWC_mean = np.array([0.010, -0.003, -0.002, -0.005]) #, np.nan, -0.002,  0.007])
#         self.DWC_std  = np.array([0.007,  0.004,  0.005,  0.004]) #, np.nan,  0.004,  0.006])
#         z_lotus  = np.array([-5.0,   -10.0,  -15.0,  -25.0]) #,  -50.0,  -75.0, -100.0])
#         # -------------------------------------

#         self.uobs =  self.DWC_mean[np.newaxis, :]
#         self.vobs = -self.CWC_mean[np.newaxis, :]
#         self.nobs =  self.CWC_mean.shape[0]
#         self.zobs = z_lotus

    

#     def forecast_score(self, 
#                        file_name: str, 
#                        ne: int = 100):
        
#         # --- Get var data ---
#         z, t, um, uv = self.get_var(file_name, 'ue', cplx=True)
        
#         # --- Interpolation at z_lotus ---
#         um_int = scp.interpolate.interp1d(z, um)(self.zobs)
#         uv_int = scp.interpolate.interp1d(z, uv)(self.zobs)
        
#         # --- MSE and VAR ---
#         MSE = (um_int.real - self.uobs)**2 + (um_int.imag - self.vobs)**2
#         VAR = uv_int.real + uv_int.imag

#         MSE = np.sum(MSE, axis = 1) / (2 * self.nobs)
#         VAR = np.sum(VAR, axis = 1) / (2 * self.nobs)

#         # --- Forecast score ---
#         cff = (ne + 1) / (ne-1)
#         EPS = np.abs(MSE - cff * VAR)
#         EPS /= EPS[0]

#         # --- Filtering ---
#         t_days = t/86400
#         dt = t_days[1] - t_days[0] # days
#         fs = 1./dt                 # samp. freq.
#         fc = 1.                    # cutoff freq. (1/days)
#         wid = fc/(fs/2)
#         num_filt, denom_filt = scp.signal.butter(5, wid, 'low')

#         return t_days, scp.signal.filtfilt(num_filt, denom_filt, EPS)
            
#     def plot_mixing(self, 
#                       files_dict: dict, 
#                       hmax: float = 50.0, 
#                       figsize: tuple = (4, 5),
#                       colors = None):
        
#         fig, ax = plt.subplots(1, 1, figsize = figsize, sharey = True)

#         kplot = 0
#         for label, file in files_dict.items():
#             path = self.output_folder + file
#             z, azm_m, azm_v = self.get_var(path, 'az_m', kt = -1, cplx = False)
            
#             if colors is not None:
#                 color = colors[kplot]
#             else:
#                 color = f'C{kplot}'
            
#             # u-profile
#             ax.plot(azm_m, z, label = label, color = color)
#             ax.fill_betweenx(z, azm_m - azm_v**0.5,
#                                 azm_m + azm_v**0.5,
#                                 alpha = 0.1, color = color)
            
#             kplot += 1
            
#         ax.set_title(r'$a_{zz}^m(z)$')
#         ax.set_xlabel(r'$[m^2.s^{-1}]$')
#         ax.set_ylabel('z (m)')
#         ax.set_ylim(-hmax, 1)
#         ax.grid()

#         handles, labels = ax.get_legend_handles_labels()
#         fig.legend(handles, labels, loc='center left', bbox_to_anchor=(1.0, 0.5), ncol=1)
#         fig.tight_layout()
        
        
# if __name__ == "__main__":
#     NENS = 1000
#     lotus = LotusTools("/home/mnex/Desktop/1DV_FD/experiments/2_Lotus_stationary/output/")

#     files_closure_impact = {r'$k-kl$' : f'k-kl_ne{NENS}.nc',}
#     lotus.plot_profiles(files_closure_impact, figsize=(8,4))
#     plt.show()