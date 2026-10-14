import torch
from taia1dvfd.core import Model, pcr

def tracer_equation_step(model: Model,
                         tracer:     torch.Tensor,
                         tracer_ref: float,
                         stochastic: bool,
                         bc_top:     str,
                         bc_btm:     str,
                         value_top:  torch.Tensor | None = None,
                         value_btm:  torch.Tensor | None = None,
                         flux_top:   float = 0.0,
                         flux_btm:   float = 0.0,
                         ) -> torch.Tensor:
    """
    One linear-implicit Euler-Maruyama step of an LU-transported tracer on cell centres,
    in banded form (no dense (ne, nz, nz) matrices).

      d C = dz( a_t dz C ) dt  -  sigma_t dz C . dB      (face-based noise, see below)

    Solved for the anomaly C - tracer_ref (float32 roundoff scales with the anomaly).

      tracer     : (ne, nz, 1) current values (not modified)
      bc_top/btm : "neumann" (flux) or "dirichlet" (value imposed in the boundary cell)
      value_*    : (ne,) or (ne, 1) Dirichlet values, in tracer units (not anomalies)
      flux_*     : kinematic fluxes w'C' at the boundary faces, positive upward
    Returns the new tracer (ne, nz, 1).
    """
    state = model.state
    cfg   = model.config
    grd   = model.grid

    # --- Variables (ne, nz) on centres, (ne, nz+1) on faces ---
    azt   = state.az_t[..., 0]
    sig_t = state.sig_z_t[..., 0]
    dB_w  = state.dB_w[..., 0]

    C_old = tracer[..., 0] - tracer_ref

    # --- Parameters ---
    dz = grd.dz
    dt = cfg.dt

    # ==================================================================
    # --- Diffusion operator D = dz( a_t dz . ), zero flux at both ends ---
    # ==================================================================
    a_face = azt[:, 1:-1]                             # (ne, nz-1): interior faces 1..nz-1

    D_up = a_face / dz**2                             # D[i, i+1]
    D_lo = a_face / dz**2                             # D[i+1, i]
    D_dg = torch.zeros_like(C_old)
    D_dg[:, :-1] -= D_up
    D_dg[:, 1:]  -= D_lo

    # --- L = I - dt D ---
    L_lo = -dt * D_lo
    L_up = -dt * D_up
    L_dg = 1.0 - dt * D_dg

    # ==================================================================
    # --- RHS ---
    # ==================================================================
    rhs = C_old.clone()

    # --- Stochastic transport (face-based, discrete fluctuation-dissipation balance) ---
    if stochastic:
        dCdz_w = torch.zeros_like(azt)                          # 0 at the boundary faces
        dCdz_w[:, 1:-1] = (C_old[:, :-1] - C_old[:, 1:]) / dz   # true dC/dz at interior faces
        noise_w = sig_t * dCdz_w * dB_w
        rhs = rhs - (noise_w[:, :-1] + noise_w[:, 1:]) / 2**0.5

    # ==================================================================
    # --- Boundary conditions ---
    # ==================================================================
    if bc_top == "neumann":
        rhs[:, 0] += -dt * flux_top / dz
    elif bc_top == "dirichlet":
        L_dg[:, 0] = 1.0
        L_up[:, 0] = 0.0
        rhs[:, 0]  = value_top.reshape(-1) - tracer_ref

    if bc_btm == "neumann":
        rhs[:, -1] += dt * flux_btm / dz
    elif bc_btm == "dirichlet":
        L_dg[:, -1] = 1.0
        L_lo[:, -1] = 0.0
        rhs[:, -1]  = value_btm.reshape(-1) - tracer_ref

    # ==================================================================
    # --- Solve ---
    # ==================================================================
    C_new = pcr(L_lo, L_dg, L_up, rhs.unsqueeze(-1))
    return C_new + tracer_ref
