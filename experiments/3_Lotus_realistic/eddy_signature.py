import numpy as np
import matplotlib

SAVE = True
matplotlib.use("Agg" if SAVE else "TkAgg")

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from mpl_toolkits.mplot3d.art3d import Line3DCollection

import sys
sys.path.append("/home/mnex/Desktop/1DV_FD/LOTUS")
from mld_obs import met_start_day, read_lotus_vel

MET_FILE = '/home/mnex/Desktop/1DV_FD/LOTUS/obs/lotusmet.epic'
VEL_FILE = '/home/mnex/Desktop/1DV_FD/LOTUS/obs/lotusnsvel.epic'
N_DAYS   = 160
STRIDE   = 4      # 15-min data -> 1 frame per hour
FPS      = 30

# --- Observations ---
t0_day = met_start_day(MET_FILE)
t_vel, z_vel, u_obs, v_obs = read_lotus_vel(VEL_FILE, t0_day, N_DAYS)
t_vel, u_obs, v_obs = t_vel[::STRIDE], u_obs[::STRIDE], v_obs[::STRIDE]
nt, nz = u_obs.shape

# --- Frame-independent quantities (computed once) ---
z      = -np.asarray(z_vel, dtype=float)
zmin   = z.min()
zbot   = np.full(nz, zmin)
vmax   = 1.1 * np.nanmax(np.hypot(u_obs, v_obs))
cm     = plt.get_cmap('cividis')
norm   = plt.Normalize(zmin, 0.0)
colors = cm(norm(z))

HEAD_RATIO = 0.15                 # head length / arrow length
HEAD_ANGLE = np.deg2rad(25)
c, s = np.cos(HEAD_ANGLE), np.sin(HEAD_ANGLE)
origin = np.column_stack([np.zeros(nz), np.zeros(nz), z])

def arrow_segments(u, v):
    """Shaft + 2 head barbs per sensor, all horizontal -> (3*nz, 2, 3)."""
    tip  = np.column_stack([np.nan_to_num(u), np.nan_to_num(v)])  # missing sensor -> zero-length
    back = -HEAD_RATIO * tip
    b1 = tip + np.column_stack([c*back[:, 0] - s*back[:, 1], s*back[:, 0] + c*back[:, 1]])
    b2 = tip + np.column_stack([c*back[:, 0] + s*back[:, 1], -s*back[:, 0] + c*back[:, 1]])
    tip3 = np.column_stack([tip, z])
    return np.concatenate([
        np.stack([origin, tip3], axis=1),
        np.stack([tip3, np.column_stack([b1, z])], axis=1),
        np.stack([tip3, np.column_stack([b2, z])], axis=1),
    ])

# --- Figure, built once ---
fig, ax = plt.subplots(figsize=(6, 7), subplot_kw={'projection': '3d'})

ax.plot([0, 0], [0, 0], [zmin, 0], color='0.4', lw=0.8)          # mooring line

arrows = Line3DCollection(arrow_segments(u_obs[0], v_obs[0]), colors=np.tile(colors, (3, 1)), linewidths=1.6)
ax.add_collection3d(arrows)
ax.set_xlim(-vmax, vmax)
ax.set_ylim(-vmax, vmax)
ax.set_zlim(zmin, 0)
ax.set_box_aspect((1, 1, 1.6))
ax.view_init(elev=20, azim=-60)
ax.set_xlabel(r'$u$ (east) [m/s]')
ax.set_ylabel(r'$v$ (north) [m/s]')
ax.set_zlabel('Depth [m]')
title = ax.set_title('')

def update(kt):
    u, v = u_obs[kt], v_obs[kt]
    arrows.set_segments(arrow_segments(u, v))
    title.set_text(f'LOTUS — day {t_vel[kt]:.2f}')
    return arrows, title # hodo, dots, title

anim = FuncAnimation(fig, update, frames=nt, interval=1000 / FPS, blit=False)

if SAVE:
    anim.save('lotus_vel.mp4', writer=FFMpegWriter(fps=FPS, bitrate=3000), dpi=150,
              progress_callback=lambda i, n: print(f'\rframe {i+1}/{n}', end=''))
    print()
else:
    plt.show()