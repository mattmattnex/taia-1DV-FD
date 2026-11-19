import torch
import numpy as np

from taia1dvfd.core  import Config, Grid, State, Model
from taia1dvfd.utils import SaveNetCDF, StepProfiler, Animator, Panel

import taia1dvfd.physics as phy
import taia1dvfd.mixing  as mix
import taia1dvfd.equation_of_state as eos
import taia1dvfd.biogeochem as bio

cfg = Config()

cfg.ne   = 100
cfg.nz   = 256
cfg.Hmin = 0.0
cfg.Hmax = 300.0
cfg.lat  = 34.98
cfg.device = "cuda:0"

# --- Time ---
cfg.dt   = 600
cfg.tend = 3600.0 * 24 * 30

# --- Wind ---
cfg.wind_mean = 5.0 + 0.0j
cfg.wind_std  = [0.0, 1.0]

# --- Wave ---
cfg.wave_mixing = True
cfg.tke_wave_mix = False
cfg.gls_wave_mix = False

# --- Stochastic bool ---
cfg.temp_stochastic = False
cfg.salt_stochastic = False
cfg.vel_stochastic  = True

# --- EOS ---
cfg.Tsurf  = 20
cfg.Tbotm  = 20
cfg.Ssurf = 38
cfg.Sbotm = 38

# --- TKE model ---
cfg.tke_pdl = True        
cfg.tke_ebb = 67.83
cfg.tke_azm_min = 1.2e-4
cfg.tke_azt_min = 1.2e-5

# --- BGC coupling ---
cfg.nbgc = 15

# --- Finalization ---
cfg.save_folder = "/home/mnex/Desktop/1DV_FD/experiments/4_Eco3M_DYFAMED/output/"
cfg.file_name = f"eco3M_ne{cfg.ne}.nc"
cfg.var_to_save = ['ue', 'mld', 'N2', 'temp']
cfg.finalize()

# --- Model objects ---
grid     = Grid(cfg)
state    = State(cfg)
model    = Model(cfg, grid, state)
savefile = SaveNetCDF(model)
prof     = StepProfiler()

# --- Initializations ---
phy.noise_LU_init(model)
phy.wind_ornstein_init(model)

mix.tke_init(model)
mix.mld_init(model)

eos.temperature_equation_init(model)
eos.salinity_equation_init(model)
eos.eos_init(model)

# --- Initialize Eco3M ---
Eco3M = bio.bgc_coupling_init(model)

DAY = 86400.0
def net_shortwave(t, qmax=800.0, albedo=0.06, daylight=0.5):
    phase = (t % DAY) / DAY                     # 0 at midnight, 0.5 at noon
    arg = (phase - 0.5) / (0.5 * daylight)
    if torch.abs(arg) >= 1.0:
        return 0.0
    return (1.0 - albedo) * qmax * torch.cos(arg * 0.5 * np.pi)

savefile.init_file(model.state)
savefile.save(model.state)

# --- Animation ---
anim = Animator(model, every=10, panels=[
Panel([("ue.real", r"$u_e$"), ("ue.imag", r"$v_e$")],
        title="Ekman velocities",
        xlabel=r"[m s$^{-1}$]",
        xlim=(-0.1, 0.1)),

Panel([("az_m", r"$a_z^m$"), ("az_t", r"$a_z^t$")],
        title="Viscosity / diffusivity",
        xlabel=r"[m$^2$ s$^{-1}$]"),

Panel([(bio.bgc_var("mid_DIN_N",   1e6), r"DIN"),
        (bio.bgc_var("mod_MODLN_N", 1e6), r"DON"),
        (bio.bgc_var("mop_MOP_N",   1e6), r"PON")],
        title="Nitrogen pools",
        xlabel=r"[$\mu$mol N L$^{-1}$]",
        xscale="log"),

Panel([(bio.bgc_var("phy_PHY_C",   1e6), r"phyto C"),
        (bio.bgc_var("bac_BAC_C",   1e6), r"bact C"),
        (bio.bgc_var("zoo_COP_C",   1e6), r"zoo C")],
        title="Living carbon",
        xlabel=r"[$\mu$mol C L$^{-1}$]",
        xscale="log"),

Panel([(bio.bgc_var("phy_PHY_Chl", 1e6), r"Chl $a$"),
        (bio.eparz(Eco3M), r"PAR")],
        title="Chlorophyll / light",
        xlabel=r"[$\mu$g L$^{-1}$]  ,  [W m$^{-2}$]",
        xscale="log"),
])

# ================== Main Loop ==================
for kt in range(cfg.nt):
    prof.start()
    
    # --- Wind ---
    phy.wind_ornstein_step(model)
    u_star = phy.compute_ustar(model).squeeze(-1)
    prof.mark("Wind")

    # --- Vertical Mixing ---
    mix.tke_step(u_star, model)
    prof.mark("Mix")

    # --- Tracers and EOS ---
    eos.temperature_equation_step(model)
    eos.salinity_equation_step(model)
    eos.eos_step(model)
    prof.mark("T-S")

    # --- Momentum ---
    phy.noise_LU_step(model)
    phy.momentum_equation_step(model)
    prof.mark("Mom")
    
    # --- Eco3M ---
    irrad = net_shortwave(model.state.time)
    bio.bgc_coupling_step(model, Eco3M, irrad)
    prof.mark("Bio")

    # --- End step ---
    model.state.time += model.config.dt
    savefile.save(model.state)
    prof.mark("Save")
    
    if kt%5 == 0:
        print(prof.format_step(kt, cfg.nt))
    anim.update(model, kt)

print(prof.format_summary())
savefile.close()
anim.close()