import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import scipy as scp

class KatoPhillipsAnalyser:
    def __init__(self, folder: str = "/home/mnex/Desktop/1DV_FD/experiments/1_Kato_Phillips/output/"):
        self.folder = folder
        self.labels = {
            'k-kl':      r'$k\text{-}kl$',
            'k-epsilon': r'$k\text{-}\epsilon$',
            'k-omega':   r'$k\text{-}\omega$',
            'gen':       r'$gen$',
            'tke':       r'$TKE$'
        }
        
    def plot_MLD(self,
        closures = ('k-kl', 'k-epsilon', 'k-omega', 'gen', 'tke'),
        nens = 500,
        dts = (36, 360, 3600),
        nzs = (512, 256),
        title = None):
        
        # Accept scalars as well as sequences
        closures = self._as_list(closures)
        dts = self._as_list(dts)
        nzs = self._as_list(nzs)

        fig, axes = plt.subplots(
            len(nzs), len(dts),
            figsize=(int(len(dts)*3 + 2), int(len(nzs)*2 + 1)),
            sharex=True, sharey=True,
            squeeze=False,          # <- always a 2D array, even for 1xN, Nx1, 1x1
        )
        if title is not None:
            fig.suptitle(title)

        for i, nz in enumerate(nzs):
            for j, dt in enumerate(dts):
                ax = axes[i, j]

                for key in closures:
                    cut = 6 if dt == 36 else 2
                    t, mld_m, mld_v = self._extract_from_file(f'{key}_dt{dt}_nz{nz}_ne{nens}.nc', cut=cut)

                    std = mld_v ** 0.5
                    line, = ax.plot(t, -mld_m, label=self.labels[key])
                    ax.fill_between(t, -(mld_m - std), -(mld_m + std),
                                    alpha=0.2, color=line.get_color())

                ax.plot(t, self._ref_MLD(3600*t), 'k--', label=r'$MLD(t)$')
                ax.set_ylim(-35, 0)
                ax.grid(True, alpha=0.4)
                if j == 0:
                    ax.set_ylabel(r'$z$ [$m$]')
                if i == len(nzs) - 1:
                    ax.set_xlabel(r'$t$ [$h$]')

        for j, dt in enumerate(dts):
            axes[0, j].annotate(rf'$\Delta t = {dt}\,$s', xy=(0.5, 1.1),
                                xycoords='axes fraction',
                                va='center', ha='center', fontsize=12,
                                bbox=dict(boxstyle='round', facecolor='peachpuff',
                                        edgecolor='0.7'))

        for i, nz in enumerate(nzs):
            axes[i, -1].annotate(rf'$n_z = {nz}$', xy=(1.045, 0.5),
                                xycoords='axes fraction',
                                va='center', ha='left', fontsize=12,
                                bbox=dict(boxstyle='round', facecolor='lemonchiffon',
                                        edgecolor='0.7'))

        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='lower center', ncol=len(labels),
                bbox_to_anchor=(0.5, -0.01))
        fig.tight_layout(rect=[0, 0.04, 1, 1])
        return fig, axes;
    
    def plot_forecast(self,
        closures = ('k-kl', 'k-epsilon', 'k-omega', 'gen', 'tke'),
        nens = 500,
        dts = (36, 360, 3600),
        nzs = (512, 256),
        title = r'Forecast score $\varepsilon$'):

        return self._plot_grid(
            compute=lambda t, m, v: self._forecast_score(t, m, v, nens),
            ylabel=r'$\varepsilon$',
            title=title, closures=closures, nens=nens, dts=dts, nzs=nzs)

    def plot_spread(self,
        closures = ('k-kl', 'k-epsilon', 'k-omega', 'gen', 'tke'),
        nens = 500,
        dts = (36, 360, 3600),
        nzs = (512, 256),
        title = r'Ensemble spread $\sigma$'):

        return self._plot_grid(
            compute=lambda t, m, v: (t, v ** 0.5),
            ylabel=r'$\sigma$',
            title=title, closures=closures, nens=nens, dts=dts, nzs=nzs)

    
    # Helpers
    def _plot_grid(self, compute, ylabel, title, closures, nens, dts, nzs):
        """Grid of subplots (rows: nz, cols: dt), one curve per closure.
        `compute(t, mld_m, mld_v) -> (x, y)` defines what is plotted."""
        
        closures = self._as_list(closures)
        dts = self._as_list(dts)
        nzs = self._as_list(nzs)

        fig, axes = plt.subplots(
            len(nzs), len(dts),
            figsize=(int(len(dts)*3 + 2), int(len(nzs)*2 + 1)),
            sharex=True, squeeze=False,
        )
        if title is not None:
            fig.suptitle(title)

        for i, nz in enumerate(nzs):
            for j, dt in enumerate(dts):
                ax = axes[i, j]

                for key in closures:
                    cut = 6 if dt == 36 else 2
                    t, mld_m, mld_v = self._extract_from_file(
                        f'{key}_dt{dt}_nz{nz}_ne{nens}.nc', cut=cut)
                    x, y = compute(t, mld_m, mld_v)
                    ax.plot(x, y, label=self.labels[key])

                ax.grid(True, alpha=0.4)
                if j == 0:
                    ax.set_ylabel(ylabel)
                if i == len(nzs) - 1:
                    ax.set_xlabel(r'$t$ [$h$]')

        for j, dt in enumerate(dts):
            axes[0, j].annotate(rf'$\Delta t = {dt}\,$s', xy=(0.5, 1.1),
                                xycoords='axes fraction',
                                va='center', ha='center', fontsize=12,
                                bbox=dict(boxstyle='round', facecolor='peachpuff',
                                          edgecolor='0.7'))

        for i, nz in enumerate(nzs):
            axes[i, -1].annotate(rf'$n_z = {nz}$', xy=(1.045, 0.5),
                                 xycoords='axes fraction',
                                 va='center', ha='left', fontsize=12,
                                 bbox=dict(boxstyle='round', facecolor='lemonchiffon',
                                           edgecolor='0.7'))

        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='lower center', ncol=len(labels),
                   bbox_to_anchor=(0.5, -0.01))
        fig.tight_layout(rect=[0, 0.04, 1, 1])
        return fig, axes
    
    def _extract_from_file(self, file, cut = 0):
        with nc.Dataset(self.folder + file) as data:
            t     = data['time'][cut:] / 3600
            mld_m = data['mld_mean'][cut:]
            mld_v = data['mld_var'][cut:]
        return t, mld_m, mld_v
    
    def _ref_MLD(self, t):
        N2 = np.sqrt(1e-4)
        ust = 0.01
        return -1.05 * ust/np.sqrt(N2) * np.sqrt(t)
    
    def _forecast_score(self, time, var_m, var_v, nr):
        # Ref MLD
        ref = self._ref_MLD(time * 3600)
        
        # Forecast score    
        MSE = (var_m - ref)**2
        VAR = (nr + 1) / (nr - 1) * var_v
        EPS = np.abs(MSE-VAR)
        
        # Coeffs for filtering
        dt = time[1] - time[0] # hours
        fs = 5./(dt) # samp. freq.
        fc = 1. # cutoff freq. (1/hours)
        wid = fc/(fs/2)
        num_filt, denom_filt = scp.signal.butter(5, wid, 'low') # 5th order Butterworth
        
        EPS = scp.signal.filtfilt(num_filt, denom_filt, EPS)
        return time, EPS
    
    @staticmethod
    def _as_list(x):
        """Return x as a list; a bare scalar or string becomes a one-element list."""
        if isinstance(x, str) or np.isscalar(x):
            return [x]
        return list(x)