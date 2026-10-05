import torch
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter

from taia1dvfd.core  import Config, Grid, State, Model, Operators
from taia1dvfd.utils import SaveNetCDF, StepProfiler

import taia1dvfd.physics as phy
import taia1dvfd.mixing  as mix
import taia1dvfd.equation_of_state as eos

cfg = Config()

cfg.ne   = 100
cfg.nz   = 256
cfg.Hmin = 1.0
cfg.Hmax = 50.0
cfg.lat  = 0.0
cfg.device = "cuda:0"

# --- Time ---
cfg.dt   = 36.0
cfg.tend = 3600.0 * 30

# --- Stochastic bool ---
cfg.temp_stochastic = False
cfg.salt_stochastic = False
cfg.vel_stochastic  = True

# --- EOS ---
cfg.Tsurf  = 13
cfg.Tbotm  = cfg.Tsurf - (1e-4 / (cfg.g * cfg.beta_T)) * (cfg.Hmax - cfg.Hmin)
cfg.Ssurf = 38
cfg.Sbotm = 38

# --- TKE model ---
cfg.tke_azm_min = 1.2e-4
cfg.tke_azt_min = 1.2e-5
cfg.tke_ebb = 67.83 # 3.75 else no wave
cfg.tke_pdl  = True

# --- Finalization ---
cfg.finalize()

# --- Model objects ---
grid   = Grid(cfg)
state  = State(cfg)
operators = Operators(cfg, grid)
model    = Model(cfg, grid, state, operators)
prof    = StepProfiler()

# --- Initializations ---
phy.noise_LU_init(model)

mix.tke_init(model)
mix.mld_init(model)

eos.temperature_equation_init(model)
eos.salinity_equation_init(model)
eos.eos_init(model)

# --- Animation setup ---
anim_every = 10
frames, frames_t = [], []

# ================== Main Loop ==================
for kt in range(cfg.nt):
    prof.start()
    
    # --- Wind forcing ---
    model.state.ua[:] = 0.01 / cfg.sqrt_cdrho
    u_star = phy.compute_ustar(model).squeeze(-1)
    prof.mark("Wind")
    
        # --- Vertical Mixing ---
    mix.tke_step(u_star, model)
    prof.mark("Mix")
    
        # --- Update mixing operators ---
    operators.update_momentum_diffusion(state.az_m)
    operators.update_tracers_diffusion(state.az_t)
    prof.mark("Ops")
    
    # --- EOS ---
    eos.temperature_equation_step(model)
    eos.salinity_equation_step(model)
    eos.eos_step(model)
    prof.mark("eos")

    # --- MLD ---
    mix.mld_step(model)
    prof.mark("mld")

    # --- Physics and stochastic
    phy.noise_LU_step(model)        
    phy.momentum_equation_step(model)
    prof.mark("Phy")
    
    # --- Update time ---
    model.state.time += model.config.dt
    
    # --- Store animation frame ---
    if kt % anim_every == 0:
        with torch.no_grad():
            ue = state.ue.squeeze(-1)
            stats = torch.stack([ue.real.mean(0), ue.real.std(0), ue.imag.mean(0), ue.imag.std(0)])
        frames.append(stats.cpu().numpy())
        frames_t.append(float(state.time))
    
    if kt%5 == 0:
        print(prof.format_step(kt, cfg.nt))

# ================== MP4 animation ==================
frames = np.stack(frames)                    # (nframes, 4, nz)
nz = frames.shape[-1]
z_src = grid.z_c if grid.z_c.shape[1] == nz else grid.z_w
z = z_src[0, :, 0].cpu().numpy()

def envelope(m, s):
    return np.column_stack([np.concatenate([m - s, (m + s)[::-1]]), np.concatenate([z, z[::-1]])])

fig, ax = plt.subplots(figsize=(3, 5))
ax.set_xlabel(r'[$m.s^{-1}$]')
ax.set_ylabel('Depth (m)')
ax.set_title('Elapsed : 0.0 h')
ax.grid(True)

# Fixed x-limits covering the whole run (mean ± std)
lim = 1.1 * np.max(np.abs(frames[:, 0::2]) + frames[:, 1::2])
ax.set_xlim(-0.01, lim)
ax.set_ylim(z.min(), z.max())

um, us, vm, vs = frames[0]
line_u, = ax.plot(um, z, label=r'$u_e$')
fill_u = ax.fill_betweenx(z, um - us, um + us, alpha=0.3, color=line_u.get_color())
ax.legend(loc='lower right')
fig.tight_layout()

writer = FFMpegWriter(fps=30, bitrate=3000)
with writer.saving(fig, "/home/mnex/Desktop/1DV_FD/experiments/1_Kato_Phillips/ekman_velocities.mp4", dpi=150):
    for i, (um, us, vm, vs) in enumerate(frames):
        line_u.set_xdata(um)
        fill_u.set_verts([envelope(um, us)])
        ax.set_title(f'Elapsed : {frames_t[i] / 3600:.2f} h')
        writer.grab_frame()

plt.close(fig)
print(f"Saved ekman_velocities.mp4 ({len(frames)} frames)")