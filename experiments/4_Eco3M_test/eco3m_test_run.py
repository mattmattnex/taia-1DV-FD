import torch
import numpy as np

from taia1dvfd.core  import Config, Grid, State, Model
from taia1dvfd.utils import SaveNetCDF, StepProfiler, Animator, Panel

import taia1dvfd.physics as phy
import taia1dvfd.mixing  as mix
import taia1dvfd.equation_of_state as eos
import taia1dvfd.biogeochem as bio

ANIM = False

# ======================================================================
cfg = Config()

cfg.ne = 20
cfg.nz = 64
cfg.Hmin = 0.0
cfg.Hmax = 150.0
cfg.lat = 34.98
cfg.device = "cuda:0"

cfg.dt = 1800.0
cfg.tend = 3600.0 * 24 * 720

# Moderate, steady wind
cfg.wind_mean = 4.0 + 0.0j
cfg.wind_std = [0.5, 0.5]

cfg.wave_mixing = True
cfg.tke_wave_mix = False
cfg.gls_wave_mix = False

# On for the ensemble
cfg.temp_stochastic = True
cfg.salt_stochastic = True
cfg.vel_stochastic  = True
cfg.bgc_stochastic  = False

cfg.Tsurf = 14.0
cfg.Tbotm = 10.0
cfg.Ssurf = 38.0
cfg.Sbotm = 38.0

cfg.tke_pdl = True
cfg.tke_ebb = 67.83
cfg.tke_azm_min = 1.2e-4
cfg.tke_azt_min = 1.2e-5

cfg.nbgc = 15

cfg.save_folder = "/home/mnex/Desktop/1DV_FD/experiments/4_Eco3M_test/output/"
cfg.file_name = "eco3M_minimal.nc"
cfg.var_to_save = [
    'mld', 'N2', 'temp',
    'bgc.phy_PHY_C', 'bgc.phy_PHY_Chl', 'bgc.phy_PHY_cell',
    'bgc.mid_DIN_N', 'bgc.mod_MODLN_N', 'bgc.mop_MOP_N',
    'bgc.bac_BAC_C', 'bgc.bac_BAC_N', 'bgc.zoo_COP_C', 'bgc.zoo_COP_N',
    'bgc.phy_PHY_N',
]
cfg.finalize()

# ======================================================================
grid = Grid(cfg)
state = State(cfg)
model = Model(cfg, grid, state)
savefile = SaveNetCDF(model)
prof = StepProfiler(sync=True)

phy.noise_LU_init(model)
phy.wind_ornstein_init(model)
mix.tke_init(model)
mix.mld_init(model)
eos.temperature_equation_init(model)
eos.salinity_equation_init(model)
eos.eos_init(model)

Eco3M = bio.bgc_coupling_init(model)

# ======================================================================
# Forcing: diel cycle modulated by a seasonal envelope.
# Day 0 is 1 January, solstice at day 172
# ======================================================================
DAY = 86400.0
YEAR = 365.0 * DAY
Q_NS = torch.full((cfg.ne,), -164.0, **cfg.farr_kwargs)

def net_shortwave(t, albedo=0.06, daylight=0.5):
    """Net surface shortwave in W/m2"""
    t = float(t)
    qmax = 550.0 + 350.0 * np.cos(2.0 * np.pi * (t - 172.0 * DAY) / YEAR)

    phase = (t % DAY) / DAY
    arg = (phase - 0.5) / (0.5 * daylight)
    if abs(arg) >= 1.0:
        return 0.0
    
    return (1.0 - albedo) * qmax * np.cos(arg * 0.5 * np.pi)


# ======================================================================
# Invariants to monitor
# ======================================================================
N_VARS = ("zoo_COP_N", "phy_PHY_N", "bac_BAC_N",
          "mop_MOP_N", "mod_MODLN_N", "mid_DIN_N")
iN = [Eco3M.index[v] for v in N_VARS]

Q_MIN_PHY_C, Q_MAX_PHY_C = 1.26e-13, 3.78e-13
i_phyC = Eco3M.index["phy_PHY_C"]
i_cell = Eco3M.index["phy_PHY_cell"]

dz = float(grid.dz)


def invariants():
    b = state.bgc
    n_tot = (b[:, :, iN].sum(dim=(1, 2)) * dz).mean().item()
    q = b[:, :, i_phyC] / b[:, :, i_cell]
    qn = ((q - Q_MIN_PHY_C) / (Q_MAX_PHY_C - Q_MIN_PHY_C))
    return n_tot, qn.min().item(), qn.max().item(), b.min().item()


N0, _, _, _ = invariants()
print(f"initial column N = {N0:.6e} mol N m-2\n")


savefile.init_file(model.state)
savefile.save(model.state)

# ======================================================================
# Animation
# ======================================================================
if ANIM:
    anim = Animator(model, every=10, panels=[
    Panel([("ue.real", r"$u_e$"), 
        ("ue.imag", r"$v_e$")],
            title="Ekman velocities",
            xlabel=r"[m s$^{-1}$]",
            xlim=(-0.1, 0.1)),

    Panel([("az_m", r"$a_z^m$"), 
        ("az_t", r"$a_z^t$")],
            title="Viscosity / diffusivity",
            xlabel=r"[m$^2$ s$^{-1}$]"),

    Panel([(bio.bgc_var("mid_DIN_N",   1e6), r"DIN"),
        (bio.bgc_var("mod_MODLN_N", 1e6), r"DON"),
        (bio.bgc_var("mop_MOP_N",   1e6), r"PON")],
            title="Nitrogen pools",
            xlabel=r"[$\mu$mol N L$^{-1}$]"),

    Panel([(bio.bgc_var("phy_PHY_C",   1e6), r"phyto C"),
        (bio.bgc_var("bac_BAC_C",   1e6), r"bact C"),
        (bio.bgc_var("zoo_COP_C",   1e6), r"zoo C")],
            title="Living carbon",
            xlabel=r"[$\mu$mol C L$^{-1}$]"),

    Panel([(bio.bgc_var("phy_PHY_Chl", 1e6), r"Chl $a$")],
            title="Chlorophyll",
            xlabel=r"[$\mu$g L$^{-1}$]  ,  [W m$^{-2}$]"),
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
    
    # --- Forcing ---
    irrad = net_shortwave(model.state.time)
    irrad_t = torch.as_tensor(irrad, **cfg.farr_kwargs)

    # --- Tracers and EOS ---
    eos.temperature_equation_step(model, heat_flux_ns=Q_NS, heat_flux_sw=irrad_t)
    eos.salinity_equation_step(model)
    eos.eos_step(model)
    mix.mld_step(model)
    prof.mark("Eos")

    # --- Momentum ---
    phy.noise_LU_step(model)
    phy.momentum_equation_step(model)
    prof.mark("Mom")
    
    # --- Eco3M ---
    bio.bgc_coupling_step(model, Eco3M, irrad)
    model.state.bgc[:] = eos.tracer_batch_equation_step(model, model.state.bgc, stochastic=cfg.bgc_stochastic)
    prof.mark("Eco3M")
    
    # --- End step ---
    model.state.time += model.config.dt
    savefile.save(model.state)    
    if ANIM: anim.update_plot(model, kt=kt)
    
    # --- BGC Invariants ---
    if kt % int(DAY / cfg.dt) == 0:
        print(prof.format_step(kt, cfg.nt))
        n_tot, qlo, qhi, cmin = invariants()
        drift = abs(n_tot - N0) / N0
        day = kt * cfg.dt / DAY
        flag = ""
        if drift > 1e-12:
            flag += "  N DRIFT"
        if qlo < -1e-6 or qhi > 1 + 1e-6:
            flag += "  QUOTA OUT OF BOUNDS"
        if cmin < 0:
            flag += "  NEGATIVE CONC"
        print(f"day {day:5.1f}  N drift {drift:8.2e} quota [{qlo:6.3f}, {qhi:6.3f}]  min {cmin:9.2e}{flag}")
        
print(prof.format_summary())
savefile.close()
if ANIM: anim.close()