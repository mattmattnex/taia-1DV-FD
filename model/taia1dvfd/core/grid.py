import torch
from dataclasses import dataclass, field
from .config import Config

@dataclass(slots=True)
class Grid:
    _cfg: Config
    
    dz:       torch.Tensor = field(init=False)
    inv_dz:   torch.Tensor = field(init=False)
    inv_2dz:  torch.Tensor = field(init=False)
    inv_dz2:  torch.Tensor = field(init=False)
    inv_2dz2: torch.Tensor = field(init=False)

    z_c: torch.Tensor = field(init=False)
    z_w: torch.Tensor = field(init=False)

    def __post_init__(self):
        cfg = self._cfg

        # Parameters
        nz = cfg.nz
        fkw = cfg.farr_kwargs
        Hmin, Hmax = cfg.Hmin, cfg.Hmax
        
        dz  = torch.tensor((Hmax - Hmin) / (nz-1), **fkw)
        self.z_c = -torch.linspace(Hmin,          Hmax +     dz, nz,   **fkw)
        self.z_w = -torch.linspace(Hmin - 0.5*dz, Hmax + 1.5*dz, nz+1, **fkw)

        self.dz = torch.abs(self.z_c[1] - self.z_c[0])
        
        self.inv_dz   = 1 / dz
        self.inv_2dz  = 1 / (2*dz)
        self.inv_dz2  = 1 / dz**2
        self.inv_2dz2 = 1 / (2*dz**2)
        
        self.z_c.unsqueeze_(0).unsqueeze_(-1)
        self.z_w.unsqueeze_(0).unsqueeze_(-1)
        
    def to_c(self, phi_w: torch.Tensor):
        """
        Going from shape _w (nz+1) -> _c (nz)
        """
        c = self._cfg
        
        phi_c = torch.zeros((c.ne, c.nz, 1), device = phi_w.device, dtype = phi_w.dtype)
        phi_c[:] = 0.5 * (phi_w[:, :-1, :] + phi_w[:, 1:, :])
        return phi_c
    
    def to_w(self, phi_c: torch.Tensor, zero_bc: bool = False):
        """
        Going from shape _c (nz) -> _w (nz+1)
        """
        c = self._cfg
        
        phi_w = torch.zeros((c.ne, c.nz+1, 1), device = phi_c.device, dtype = phi_c.dtype)
        phi_w[:, 1:-1, :] = 0.5 * (phi_c[:, :-1, :] + phi_c[:, 1:, :])
        
        if not zero_bc:
            phi_w[:,  0, :] = phi_w[:,  1, :]
            phi_w[:, -1, :] = phi_w[:, -2, :]
        
        return phi_w