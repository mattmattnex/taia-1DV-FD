import torch
import numpy as np
from taia1dvfd.core import Model
from .eco3m_coupler import Eco3MCoupler

def bgc_coupling_init(model: Model) -> Eco3MCoupler:
    cfg = model.config
    grd = model.grid
    
    bgc = Eco3MCoupler(libpath="Eco3M/MAKE/libeco3m.so", 
                       workdir="Eco3M/MAKE")
    
    bgc.initialize(nz = cfg.nz)
    bgc.set_dt_bio(cfg.dt)
    bgc.set_dz(np.full(cfg.nz, grd.dz.cpu().numpy()))
    
    return bgc

def bgc_coupling_step(model: Model, bgc: Eco3MCoupler):
    cfg = model.config
    state = model.state
    
    # Fetch state
    temp = state.temp.cpu().numpy()[0, :, 0]
    salt = state.salt.cpu().numpy()[0, :, 0]
    irrad = 0

    new_conc = bgc.step(bgc.CONC, temp, salt, irrad, model.state.time)
    bgc.CONC += cfg.dt * new_conc