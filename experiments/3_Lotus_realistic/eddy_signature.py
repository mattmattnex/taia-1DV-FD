import numpy as np
from scipy.ndimage import uniform_filter1d

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

import sys
sys.path.append("/home/mnex/Desktop/1DV_FD/LOTUS")
from mld_obs import met_start_day, read_lotus_vel

MET_FILE = '/home/mnex/Desktop/1DV_FD/LOTUS/obs/lotusmet.epic'
VEL_FILE = '/home/mnex/Desktop/1DV_FD/LOTUS/obs/lotusnsvel.epic'
FPS = 30

days_to_plot = range(90, 130)
average_window = 8

zlim = (-50, 0)
vmax = 1.0

def _block_mean(a, window):
    nb = a.shape[0] // window
    a = a[: nb * window].reshape(nb, window, *a.shape[1:])
    return np.nanmean(a, axis=1)

def get_days_series(t, u, v, days, window=0):
    mask = np.isin(np.floor(t).astype(int), np.asarray(list(days)))
    ts, us, vs = t[mask], u[mask], v[mask]
    ts, us, vs = (_block_mean(a, window) for a in (ts, us, vs))
    return ts, us, vs

# --- Observations ---
t0_day = met_start_day(MET_FILE)
t, z, u, v = read_lotus_vel(VEL_FILE, t0_day, 160)
t, u, v = get_days_series(t, u, v, days_to_plot, average_window)
nt, nz = u.shape

# --- Figure, built once ---
fig, ax = plt.subplots(figsize=(6, 7), subplot_kw={'projection': '3d'})
cm = plt.get_cmap('viridis')
norm = plt.Normalize(zlim[0], 0.0)

def update(kt):
    uk, vk = u[kt], v[kt]
    
    ax.clear()
    ax.set_xlim3d(-vmax, vmax)
    ax.set_ylim3d(-vmax, vmax)
    ax.set_zlim3d(*zlim)
    ax.set_xlabel(r'$u$ [m/s]')
    ax.set_ylabel(r'$v$ [m/s]')
    ax.set_box_aspect((1, 1, 1.6))
    ax.view_init(elev=90, azim=90)
    ax.set_zlabel('Depth [m]')
    ax.plot([0, 0], [0, 0], [*zlim], color='0.4', lw=0.8) # Mooring
    ax.set_title(f'LOTUS — day {t[kt]:.2f}')

    # Arrows, coloured by depth
    for ui, vi, zi in zip(uk, vk, -z):
        ax.quiver(0, 0, zi, ui, vi, 0, 
                  color=cm(norm(zi)), 
                  arrow_length_ratio=0.15, 
                  linewidth=1.6, 
                  normalize=False)
        
anim = FuncAnimation(fig, update, frames=nt, interval=1000 / FPS, blit=False)

anim.save('lotus_vel.mp4', writer=FFMpegWriter(fps=FPS, bitrate=3000), dpi=150,
            progress_callback=lambda i, n: print(f'\rframe {i+1}/{n}', end=''))