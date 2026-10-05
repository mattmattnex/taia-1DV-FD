import torch
import numpy as np

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from taia1dvfd.core import Model

class Animator:
    def __init__(self, model: Model):
        
        cfg = model.config
        grid = model.grid
        self.z_c   = grid.z_c[0, :, 0].cpu().numpy()
        self.z_w   = grid.z_w[0, :, 0].cpu().numpy()
        self.zer_c = np.zeros_like(self.z_c)
        self.zer_w = np.zeros_like(self.z_w)
        
        self.fig, self.ax = plt.subplots(1, 5, figsize = (14, 5), sharey=True)
        
        # Ekman velocities
        self.ax[0].set_title(r'Ekman velocities $\boldsymbol{u}_e$')
        self.ax[0].grid(True)
        self.ax[0].set_xlim(-0.1, 0.1)
        self.ax[0].set_ylabel('Depth (m)')
        self.ax[0].set_xlabel(r'[$m.s^{-1}$]')
        self.ue_plot  = self._gen_line_poly(self.ax[0], label = r"$u_e$")
        self.ve_plot  = self._gen_line_poly(self.ax[0], label = r"$v_e$")
        self.ax[0].legend()
        
        # Eddy visc and diff
        self.ax[1].set_title(r'Viscosity $a_{zz}^m$ and diffusivity $a_{zz}^t$')
        self.ax[1].grid(True)
        self.ax[1].set_xlim(-0.01, 0.1)
        self.ax[1].set_xlabel(r'[$m.s^{-2}$]')
        self.azm_plot = self._gen_line_poly(self.ax[1], label = r"$a_{z,m}$")
        self.azt_plot = self._gen_line_poly(self.ax[1], label = r"$a_{z,t}$")
        self.ax[1].legend()
        
        # Temperature and Salinity
        self.ax[2].set_title(r'Temperature')
        self.ax[2].grid(True)
        self.ax[2].set_xlim(19, 22)
        self.ax[2].set_xlabel(r'[$^\circ C$]')
        self.temp_plot = self._gen_line_poly(self.ax[2], label=r"$T$")
        self.ax[2].legend()

        # N2
        self.ax[3].set_title(r'Brunt-Vaisälä freq')
        self.ax[3].grid(True)
        self.ax[3].set_xlim(-1e-3, 1e-3)
        self.N2_plot = self._gen_line_poly(self.ax[3], label=r"$N^2$")
        self.ax[3].legend()

        # Density
        self.ax[4].set_title(r'Density')
        self.ax[4].grid(True)
        self.ax[4].set_xlim(1025, 1027)
        self.ax[4].set_xlabel(r'[$kg.m^{-3}$]')
        self.rho_plot = self._gen_line_poly(self.ax[4], label=r"$\rho$")
        self.ax[4].legend()
        
    def update_plot(self, model:Model):
        state = model.state
        
        self.fig.suptitle(f'Elapsed : {state.time}s')
        
        self._update_line_poly(state.ue.real, self.ue_plot)
        self._update_line_poly(state.ue.imag, self.ve_plot)
        self._update_line_poly(state.az_m,    self.azm_plot)
        self._update_line_poly(state.az_t,    self.azt_plot)
        self._update_line_poly(state.N2,      self.N2_plot)
        self._update_line_poly(state.temp,    self.temp_plot)
        self._update_line_poly(state.rho,     self.rho_plot)
        
        plt.pause(1e-5)
        
                
    # ==========================        
    # --- Internal functions ---
    # ==========================
    
    def _update_line_poly(self, var: torch.Tensor, plot_tuple: tuple) -> None:
        var_m, var_std = self._get_stats_np(var)
        
        nz = var.shape[1]
        if nz == self.z_c.shape[0]:
            z   = self.z_c
        else:
            z = self.z_w

        # Update Line
        plot_tuple[0].set_data(var_m, z)
        
        # Update poly (variance enveloppe)
        poly_x = np.concatenate([var_m - var_std, (var_m + var_std)[::-1]])
        poly_y = np.concatenate([z, z[::-1]])
        plot_tuple[1].set_verts([np.column_stack([poly_x, poly_y])])
    
    def _gen_line_poly(self, ax, label = None, color = None) -> None:
        z   = self.z_c
        zer = self.zer_c 
            
        line = ax.plot(zer, z, label = label)
        poly = ax.fill_betweenx(z, zer, zer, alpha = 0.3)
        return (line[0], poly)
        
    @staticmethod
    def _get_stats_np(var: torch.Tensor):
        var_mean = torch.mean(var.squeeze(-1), dim = 0)
        var_std  = torch.std(var.squeeze(-1),  dim = 0)
        
        return var_mean.cpu().numpy(), var_std.cpu().numpy()