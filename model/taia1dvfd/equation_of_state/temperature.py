import torch
from taia1dvfd.core import Model
from .tracer import tracer_equation_step

def temperature_equation_init(model: Model, temp_init: torch.Tensor | None = None):
    zc = model.grid.z_c
    cfg = model.config

    if temp_init is not None:
        T0 = temp_init.clone()

    else:
        T0 = torch.zeros_like(model.state.temp)
        dT = cfg.Tsurf - cfg.Tbotm
        T0[:] = cfg.Tsurf - dT * (-zc - cfg.Hmin) / (cfg.Hmax - cfg.Hmin)

    model.state.temp[:] = T0

def temperature_equation_step(model: Model,
                              heat_flux_ns: torch.Tensor | None = None,
                              heat_flux_sw: torch.Tensor | None = None,
                              temp_surf: torch.Tensor | None = None,
                              temp_botm: torch.Tensor | None = None,
                              ):
    """
    heat_flux_ns   : (ne,) non-solar heat flux INTO the ocean (W/m2) -> surface Neumann BC
    heat_flux_sw   : scalar or (ne,) net shortwave INTO the ocean (W/m2), absorbed over depth
    temp_surf/botm : Dirichlet values, only used if cfg.temp_bc_top/btm = "dirichlet"
    """
    cfg = model.config
    rho_cp = cfg.rho_o * cfg.eos_cp

    # --- Surface flux: w'T'(top) = -Q_ns / (rho0 cp), positive upward ---
    flux_top = 0.0
    if heat_flux_ns is not None:
        flux_top = -heat_flux_ns.reshape(-1) / rho_cp

    # --- Penetrating shortwave (K/s) ---
    source = None
    if heat_flux_sw is not None:
        source = shortwave_heating(model, heat_flux_sw) / rho_cp

    model.state.temp[:] = tracer_equation_step(
        model,
        tracer     = model.state.temp,
        tracer_ref = cfg.eos_T0,
        stochastic = cfg.temp_stochastic,
        bc_top     = cfg.temp_bc_top,
        bc_btm     = cfg.temp_bc_btm,
        value_top  = temp_surf,
        value_btm  = temp_botm,
        flux_top   = flux_top,
        flux_btm   = 0.0,
        source     = source,
    )

def shortwave_heating(model: Model,
                      heat_flux_sw: torch.Tensor,
                      R: float = 0.58,
                      zeta_1: float = 0.35,
                      zeta_2: float = 23.0,
                      ) -> torch.Tensor:
    """
    Shortwave absorbed in each cell (W/m3), Paulson & Simpson (1977), Jerlov type I:

      I(z) = Q_sw [ R exp(z / zeta_1) + (1 - R) exp(z / zeta_2) ],   z <= 0

    The top face is at -Hmin, so what is absorbed above it (~55 % for Hmin = 1 m)
    goes into the first cell. What passes the bottom face (~0.5 % at 100 m) leaves.
    Returns (ne, nz) or (1, nz).
    """
    grd = model.grid
    zc  = grd.z_c.reshape(-1)
    dz  = grd.dz

    # --- Faces (top -> bottom), uniform grid ---
    z_w = torch.cat((zc[:1] + dz / 2, zc - dz / 2))

    # --- Fraction of Q_sw crossing each face ---
    sw_frac = R * torch.exp(z_w / zeta_1) + (1 - R) * torch.exp(z_w / zeta_2)
    sw_frac[0] = 1.0

    absorbed = (sw_frac[:-1] - sw_frac[1:]) / dz                 # (nz,) m-1
    return torch.as_tensor(heat_flux_sw).reshape(-1, 1) * absorbed[None, :]