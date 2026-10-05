import torch
import torch.nn.functional as F

from taia1dvfd.core import Model

def mld_init(model: Model):
    state = model.state
    state.mld[:] = 0.0
    
def mld_step(model: Model):
    state = model.state
    grd = model.grid 
    
    N2 = grd.to_c(state.N2).squeeze(-1)
    z  = grd.z_w.squeeze()

    idx = N2[:, 2:].argmax(dim=1) + 2

    model.state.mld[:, 0, 0] = torch.abs(z[idx])