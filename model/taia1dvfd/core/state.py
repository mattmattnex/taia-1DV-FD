import torch
from dataclasses import dataclass, field
from .config import Config

@dataclass(slots=True)
class State:
    _cfg: Config
    
    # Time
    time: torch.Tensor = field(init=False)

    # Wind, waves and ekman vel
    ua: torch.Tensor = field(init=False)
    us: torch.Tensor = field(init=False)
    ue: torch.Tensor = field(init=False)
    
    # Noise
    dB_w : torch.Tensor = field(init=False)
    dB_c : torch.Tensor = field(init=False)
    sig_h: torch.Tensor = field(init=False)
    sig_z_m: torch.Tensor = field(init=False)
    sig_z_t: torch.Tensor = field(init=False)

    # Vertical mixing
    az_m: torch.Tensor = field(init=False)
    az_t: torch.Tensor = field(init=False)
    
    # TKE closure
    tke:   torch.Tensor = field(init=False)
    dissl: torch.Tensor = field(init=False)

    # GLS closure
    psi:  torch.Tensor = field(init=False)
    lsc:  torch.Tensor = field(init=False)
    diss: torch.Tensor = field(init=False)
    M2:   torch.Tensor = field(init=False)

    # EOS
    temp: torch.Tensor = field(init=False)
    salt: torch.Tensor = field(init=False)
    rho:  torch.Tensor = field(init=False)
    N2:   torch.Tensor = field(init=False)
    mld:  torch.Tensor = field(init=False)
    
    def __post_init__(self):
        cfg = self._cfg
        
        ne, nz = cfg.ne, cfg.nz
        ckw = cfg.carr_kwargs
        fkw = cfg.farr_kwargs
        
        # Time
        self.time = torch.tensor(0.0, **fkw)

        # Physics
        self.ua    = torch.zeros(ne,    1, 1, **ckw)
        self.us    = torch.zeros(ne, nz+1, 1, **ckw)
        self.ue    = torch.zeros(ne, nz+1, 1, **ckw)
        
        # Sigmas
        self.sig_h   = torch.zeros(ne, nz+1, 1, **ckw)
        self.sig_z_m = torch.zeros(ne, nz+1, 1, **fkw)
        self.sig_z_t = torch.zeros(ne, nz+1, 1, **fkw)
        
        # Noise
        self.dB_w  = torch.zeros(ne, nz+1, 1, **fkw)
        self.dB_c  = torch.zeros(ne,   nz, 1, **fkw)

        # Vertical mixing
        self.az_m = torch.zeros(ne, nz+1, 1, **fkw)
        self.az_t = torch.zeros(ne, nz+1, 1, **fkw)
        
        # TKE closure
        self.tke   = torch.zeros(ne, nz+1, 1, **fkw)
        self.dissl = torch.zeros(ne, nz+1, 1, **fkw)

        # EOS
        self.temp = torch.zeros(ne,   nz, 1, **fkw)
        self.salt = torch.zeros(ne,   nz, 1, **fkw)
        self.rho  = torch.zeros(ne,   nz, 1, **fkw)
        self.N2   = torch.zeros(ne, nz+1, 1, **fkw)
        self.mld  = torch.zeros(ne,    1, 1, **fkw)
        
        # GLS closure
        self.psi  = torch.zeros(ne, nz+1, 1, **fkw)
        self.lsc  = torch.zeros(ne, nz+1, 1, **fkw)
        self.diss = torch.zeros(ne, nz+1, 1, **fkw)
        self.M2   = torch.zeros(ne, nz+1, 1, **fkw)