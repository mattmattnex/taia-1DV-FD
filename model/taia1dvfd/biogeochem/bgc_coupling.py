import torch
import numpy as np
from taia1dvfd.core import Model
from .eco3m_coupler import Eco3MCoupler
from .eco3m_cini import load_cini

CINI_DIR = "/home/mnex/Desktop/1DV_FD/Eco3M/CONFIG_ECO3M/C_INI/1DV/C_INI_DYFAMED"
DEPTH_FILE = "/home/mnex/Desktop/1DV_FD/Eco3M/OUTPUTS/sim_eco3m_1DV_depth.dat"

def bgc_coupling_init(model: Model) -> Eco3MCoupler:
    cfg = model.config
    grd = model.grid
    
    bgc = Eco3MCoupler(libpath="Eco3M/MAKE/libeco3m.so", workdir="Eco3M/MAKE")
    bgc.initialize(nz = cfg.nz, ne = cfg.ne)
    bgc.set_dt_bio(cfg.dt)
    bgc.set_dz(np.full(cfg.nz, grd.dz.cpu().numpy()))
    
    if bgc.nbvar != cfg.nbgc:
        raise ValueError(f"cfg.nbgc={cfg.nbgc} but config.ini gives nbvar={bgc.nbvar}")
    
    # Load initial config
    z = -grd.z_c.squeeze().cpu().numpy()          # positive, surface first
    conc = load_cini(CINI_DIR, DEPTH_FILE, z, bgc.varnames)
    
    # Setup initial concentrations
    model.state.bgc[:] = torch.from_numpy(conc).to(cfg.device)   # all members
    cfg.bgc_names = bgc.varnames
    
    return bgc

def bgc_coupling_step(model: Model, bgc: Eco3MCoupler, irrad: float):
    cfg = model.config
    state = model.state
    
    temp = state.temp.cpu().numpy()[:, :, 0].astype(np.float64)
    salt = state.salt.cpu().numpy()[:, :, 0].astype(np.float64)
    conc = state.bgc.cpu().numpy()

    tend = bgc.step(conc, temp, salt, irrad, float(state.time))
    state.bgc += cfg.dt * torch.from_numpy(tend).to(cfg.device)