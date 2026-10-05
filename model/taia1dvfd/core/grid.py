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
        
        # Domain [-Hmax, -Hmin], nz cells of thickness dz.
        # Faces (z_w, nz+1) bound the cells; centers (z_c, nz) are midpoints.
        dz = (Hmax - Hmin) / nz
        self.z_w = -torch.linspace(Hmin, Hmax, nz + 1, **fkw)
        self.z_c = 0.5 * (self.z_w[:-1] + self.z_w[1:])

        self.dz = torch.tensor(dz, **fkw)
        self.inv_dz   = 1 / self.dz
        self.inv_2dz  = 1 / (2 * self.dz)
        self.inv_dz2  = 1 / self.dz**2
        self.inv_2dz2 = 1 / (2 * self.dz**2)
        
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