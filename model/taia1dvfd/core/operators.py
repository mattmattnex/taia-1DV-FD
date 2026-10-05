import torch
from dataclasses import dataclass, field
from .config import Config
from .grid import Grid

@dataclass(slots=True)
class Operators:
    _cfg: Config
    _grd: Grid
    
    d_dz_c:   torch.Tensor = field(init=False)
    d_dz_w:   torch.Tensor = field(init=False)
    d_dz_c2w: torch.Tensor = field(init=False)
    d_aztdz:  torch.Tensor = field(init=False)
    d_azmdz:  torch.Tensor = field(init=False)
    
    def __post_init__(self):
        cfg = self._cfg
        ne, nz = cfg.ne, cfg.nz
        fkw = cfg.farr_kwargs
        ckw = cfg.carr_kwargs
        
        self.d_dz_c   = torch.zeros((ne, nz,   nz  ), **fkw)
        self.d_dz_w   = torch.zeros((ne, nz+1, nz+1), **ckw)
        self.d_dz_c2w = torch.zeros((ne, nz+1, nz  ), **fkw)
        self.d_aztdz  = torch.zeros((ne, nz,   nz  ), **fkw)
        self.d_azmdz  = torch.zeros((ne, nz+1, nz+1), **fkw)
        
        self._compute_d_dz()
        
    def _compute_d_dz(self):
        nz = self._cfg.nz
        cff = self._grd.inv_dz

        # Cell centers
        for i in range(1, nz-1):
            self.d_dz_c[:, i, i-1] = -0.5
            self.d_dz_c[:, i, i+1] =  0.5
            
        self.d_dz_c[:, 0, 0] = -1.5
        self.d_dz_c[:, 0, 1] = 2
        self.d_dz_c[:, 0, 2] = -0.5
        self.d_dz_c[:, -1, -1] = 1.5
        self.d_dz_c[:, -1, -2] = -2
        self.d_dz_c[:, -1, -3] = 0.5
        self.d_dz_c *= -cff
        
        # Cell faces
        for i in range(1, nz):
            self.d_dz_w[:, i, i-1] = -0.5
            self.d_dz_w[:, i, i+1] =  0.5
            
        self.d_dz_w[:, 0, 0] = -1.5
        self.d_dz_w[:, 0, 1] = 2
        self.d_dz_w[:, 0, 2] = -0.5
        self.d_dz_w[:, -1, -1] = 1.5
        self.d_dz_w[:, -1, -2] = -2
        self.d_dz_w[:, -1, -3] = 0.5
        self.d_dz_w *= -cff
        
        # Cell centers to faces
        for i in range(1, nz):
            self.d_dz_c2w[:, i, i-1] = -1.0
            self.d_dz_c2w[:, i, i  ] =  1.0
        self.d_dz_c2w *= -cff

    def update_momentum_diffusion(self, azm):
        nz, cff = self._cfg.nz, self._grd.inv_dz2
        
        a_half = 0.5 * (azm[:, :-1, 0] + azm[:, 1:, 0])      # (ne, nz): a at i+1/2
        i = torch.arange(1, nz, device=azm.device)
        
        self.d_azmdz.zero_()
        self.d_azmdz[:, i, i - 1] =  a_half[:, i - 1] * cff
        self.d_azmdz[:, i, i + 1] =  a_half[:, i]     * cff
        self.d_azmdz[:, i, i]     = -(a_half[:, i - 1] + a_half[:, i]) * cff


    def update_tracers_diffusion(self, azt: torch.Tensor):
        nz  = self._cfg.nz
        cff = self._grd.inv_dz2

        i = torch.arange(1, nz - 1, device=azt.device)

        azt_up = azt[:, i,     0]
        azt_dn = azt[:, i + 1, 0]

        self.d_aztdz.zero_()
        self.d_aztdz[:, i, i - 1] = azt_up * cff
        self.d_aztdz[:, i, i + 1] = azt_dn * cff
        self.d_aztdz[:, i, i    ] = -(azt_up + azt_dn) * cff