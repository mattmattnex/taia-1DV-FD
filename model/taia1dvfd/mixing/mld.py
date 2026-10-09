import torch
import torch.nn.functional as F

from taia1dvfd.core import Model

def mld_init(model: Model):
    state = model.state
    state.mld[:] = 0.0

def mld_step(model: Model):
    """
    Mixed layer depth = depth of the maximum of N2 (entrainment interface).

    1. N2 is moved from faces (nz+1) to cell centres (nz).
    2. N2 is smoothed with a running mean over cfg.mld_l_smooth [m] to remove
       grid-scale noise. No smoothing when dz >= mld_l_smooth (coarse grids).
    3. The peak is located with argmax, then refined to sub-cell precision
       with a parabolic fit through the peak and its two neighbours.
    """
    cfg   = model.config
    state = model.state
    grd   = model.grid

    # --- Parameters ---
    nz       = cfg.nz
    dz       = (cfg.Hmax - cfg.Hmin) / nz            # Python float: no grd.dz.item() -> no GPU sync
    l_smooth = cfg.mld_l_smooth
    n_skip   = 2                                     # cells ignored right below the surface

    # --- N2 at cell centres (no temporary zeros buffer) ---
    N2_w    = state.N2[..., 0]                       # (ne, nz+1)
    N2_c    = 0.5 * (N2_w[:, :-1] + N2_w[:, 1:])     # (ne, nz)
    depth_c = -grd.z_c.view(-1)                      # (nz) depth of centres, positive downward

    # --- Running-mean smoothing over l_smooth ---
    # avg_pool1d == conv1d with a ones/K kernel, but a dedicated kernel
    # (conv1d in float64 has no cuDNN path and falls back to a slow one)
    half_width = int(round(0.5 * l_smooth / dz))     # 0 -> no smoothing

    if half_width > 0:
        kernel_size = 2 * half_width + 1
        N2_pad  = F.pad(N2_c.unsqueeze(1), (half_width, half_width), mode="replicate")
        N2_smth = F.avg_pool1d(N2_pad, kernel_size=kernel_size, stride=1).squeeze(1)  # (ne, nz)
    else:
        N2_smth = N2_c

    # --- Peak index (searched below the first n_skip cells) ---
    idx_peak = N2_smth[:, n_skip:].argmax(dim=1) + n_skip
    idx_peak = idx_peak.clamp(1, nz - 2)             # both neighbours must exist

    # --- Values around the peak: one gather for the 3 points ---
    offsets = torch.arange(-1, 2, device=idx_peak.device)                    # (-1, 0, +1)
    N2_up, N2_peak, N2_dn = N2_smth.gather(1, idx_peak.unsqueeze(1) + offsets).unbind(1)

    # --- Parabolic refinement ---
    # Vertex of the parabola through (-1, N2_up), (0, N2_peak), (+1, N2_dn),
    # in cell units, positive = deeper. Only valid for a true maximum (curvature < 0).
    curvature   = N2_up - 2.0 * N2_peak + N2_dn
    is_max      = curvature < 0.0
    curv_safe   = torch.where(is_max, curvature, -1.0)

    shift_cells = torch.where(is_max, 0.5 * (N2_up - N2_dn) / curv_safe, 0.0)
    shift_cells = shift_cells.clamp(-0.5, 0.5)

    # --- Mixed layer depth ---
    mld = depth_c[idx_peak] + shift_cells * dz
    state.mld[:, 0, 0] = mld