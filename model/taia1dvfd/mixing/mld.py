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
    dz       = grd.dz.item()
    l_smooth = cfg.mld_l_smooth
    n_skip   = 2                                     # cells ignored right below the surface
 
    # --- N2 at cell centres ---
    N2_c    = grd.to_c(state.N2)[..., 0]             # (ne, nz)
    depth_c = -grd.z_c.view(-1)                      # (nz) depth of centres, positive downward
 
    # --- Running-mean smoothing over l_smooth ---
    half_width = int(round(0.5 * l_smooth / dz))     # 0 -> no smoothing
    
    if half_width > 0:
        kernel_size = 2 * half_width + 1
        kernel  = torch.ones(1, 1, kernel_size, dtype=N2_c.dtype, device=N2_c.device) / kernel_size
        N2_pad  = F.pad(N2_c.unsqueeze(1), (half_width, half_width), mode="replicate")
        N2_smth = F.conv1d(N2_pad, kernel).squeeze(1)  # (ne, nz)
    else:
        N2_smth = N2_c
 
    # --- Peak index (searched below the first n_skip cells) ---
    idx_peak = N2_smth[:, n_skip:].argmax(dim=1) + n_skip
    idx_peak = idx_peak.clamp(1, nz - 2)             # both neighbours must exist
 
    # --- Values around the peak ---
    N2_up   = N2_smth.gather(1, (idx_peak - 1).unsqueeze(1)).squeeze(1)   # cell above (shallower)
    N2_peak = N2_smth.gather(1,  idx_peak     .unsqueeze(1)).squeeze(1)
    N2_dn   = N2_smth.gather(1, (idx_peak + 1).unsqueeze(1)).squeeze(1)   # cell below (deeper)
 
    # --- Parabolic refinement ---
    # Vertex of the parabola through (-1, N2_up), (0, N2_peak), (+1, N2_dn),
    # in cell units, positive = deeper. Only valid for a true maximum (curvature < 0).
    curvature   = N2_up - 2.0 * N2_peak + N2_dn
    is_max      = curvature < 0.0
    curv_safe   = torch.where(is_max, curvature, -torch.ones_like(curvature))
    
    shift_cells = torch.where(is_max, 0.5 * (N2_up - N2_dn) / curv_safe, torch.zeros_like(curvature))
    shift_cells = shift_cells.clamp(-0.5, 0.5)
 
    # --- Mixed layer depth ---
    mld = depth_c[idx_peak] + shift_cells * dz
    state.mld[:, 0, 0] = mld