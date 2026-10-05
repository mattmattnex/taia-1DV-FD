import math
from dataclasses import dataclass, field
from typing import Dict, Any, List
import torch

class Config:
    # --- Model parameters ---
    ne: int = 100          # Number of ensemble members
    nz: int = 256          # Number of vertical levels (T-grid)
    Hmin: float = 1.0      # Minimum depth
    Hmax: float = 256.0    # Maximum depth
    lat: float = 34.0      # Latitude (deg N)
    device: str = "cuda:0"

    # --- Time parameters ---
    dt: float = 1800.0                # Time step (s)
    tend: float = 3600.0 * 24 * 20    # End time (s)

    # --- Simulation ---
    save_folder = "/home/mnex/Desktop/1DV_FD/experiments/"
    file_name = f"output.nc"
    var_to_save = ['temp', 'mld', 'N2']

    # --- Physics ---
    vel_stochastic: bool = True
    rho_o: float = 1000.0
    rho_a: float = 1.0
    Cd: float = 1.3e-3
    g: float = 9.81
    cp: float = 3985.0
    
    # --- Wind ---
    wind_mean: complex = 5.0 + 0j
    wind_std:  float = 5.0
    wind_Ta:   float = 86400.0

    # --- Itô-Stokes ---
    us0: float = 0.068
    hs:  float = 4.775

    # --- Waves ---
    theta_mean: float = 0.0
    theta_std:  float = math.radians(5.0)
    wave_mixing: bool = True
    
    # --- KPP ---
    kpp_c1:      float = 0.4
    kpp_c2:      float = 0.7
    kpp_zt_0:    float = 0.05
    kpp_azm_min: float = 1e-4
    kpp_azt_min: float = 1e-4
    kpp_regularization: bool = False
    
    # --- TKE closure ---
    tke_ediff:     float = 0.1        # rn_ediff : C_k, viscosity coef.  avm = Ck * l * sqrt(e)
    tke_ediss:     float = 0.7        # rn_ediss : C_eps, dissipation coef.
    tke_ebb:       float = 67.83      # rn_ebb   : surface TKE input coef.
    tke_efr:       float = 0.05       # fraction of surface TKE that penetrates
    tke_htau:      float = 10.0       # e-folding depth of TKE penetration  [m]
    tke_min:       float = 1.0e-6     # rn_emin  : minimum interior TKE   [m2/s2]
    tke_surf_min:  float = 1.0e-4     # rn_emin0 : minimum surface TKE    [m2/s2]
    tke_bshear:    float = 1.0e-20    # rn_bshear: background shear       [s-2]
    tke_azm_floor: float = 2.0e-5     # floor on avm in TKE diffusion solve
    tke_rmxl0:     float = 0.04
    tke_N2_min:    float = 1e-14
    tke_pdl:       bool  = True
    tke_wave_mix:  bool  = True
    tke_azm_min:   float = 1.2e-4
    tke_azt_min:   float = 1.2e-5

    # --- GLS ---
    gls_model: str = "k-kl"
    gls_stabf: str = "CA"
    gls_bc:    str = "dirichlet"
    gls_wave_mix: bool  = True
    gls_vonKar:   float = 0.4
    gls_galp:     float = 0.53
    gls_chk:      float = 1400 / 9.81
    gls_zos_min:  float = 1e-2
    gls_zob_min:  float = 1e-2
    gls_lsc0:     float = 0.05
    gls_cmu_min:      float = 1e-10
    gls_cmu_prim_min: float = 1e-10
    gls_M2_min:   float = 1e-20
    gls_N2_min:   float = 1e-14
    gls_lsc_min:  float = 1e-8
    gls_psi_min:  float = 1e-12
    gls_azm_min:  float = 1.2e-4
    gls_azt_min:  float = 1.2e-5
    gls_tke_min:  float = 1e-6
    gls_diss_min: float = 1e-12
    
    _eps_sf: float = 1e-20
    _lim_sf: float = 1e10
    
    gls_model_params: dict = {
    'k-kl':      {'rp':  0.0, 'rm': 1.0, 'rn':  1.0,     'beta1': 0.9,   'beta2': 0.5,   'beta3p': 1.0, 'beta3m': {'KC':  2.53, 'CA':  2.62}, 'sig_tke': 1.96, 'sig_psi': 1.96},
    'k-epsilon': {'rp':  3.0, 'rm': 1.5, 'rn': -1.0,     'beta1': 1.44,  'beta2': 1.92,  'beta3p': 1.0, 'beta3m': {'KC': -0.52, 'CA': -0.63}, 'sig_tke': 1.0,  'sig_psi': 1.3},
    'k-omega':   {'rp': -1.0, 'rm': 0.5, 'rn': -1.0,     'beta1': 0.555, 'beta2': 0.833, 'beta3p': 1.0, 'beta3m': {'KC': -0.58, 'CA': -0.64}, 'sig_tke': 2.0,  'sig_psi': 2.0},
    'gen':       {'rp':  2.0, 'rm': 1.0, 'rn': -2.0/3.0, 'beta1': 1.0,   'beta2': 1.22,  'beta3p': 1.0, 'beta3m': {'KC':  0.1,  'CA':  0.05}, 'sig_tke': 0.8,  'sig_psi': 1.07}
    }
    
    gls_stabf_params: dict = {
        "KC": {"c":  [6.00, 0.32, 0.000, 0.000, 0.000, 0.000], "cb": [3.728, 0.70, 0.70, 0.00, 0.20,   0.6102]},
        "CA": {"c":  [5.00, 0.80, 1.968, 1.136, 0.000, 0.400], "cb": [5.95,  0.60, 1.00, 0.00, 0.3333, 0.72 ]},
    }

    # --- Equation of state --- #
    eos_T0: float = 20   # Reference temperature (°C)
    eos_S0: float = 38   # Reference Salinity (PSU)
    beta_T: float = 2e-4 # Thermal expansion coefficient
    beta_S: float = 1e-3 # Haline contraction coefficient
    
    # --- Salinity model ---
    salt_stochastic: bool = False
    Ssurf:  float = 20.0 # Initial surface S
    Sbotm:  float = 20.5 # Initial bottom S
    
    # --- Temperature model ---
    temp_stochastic: bool = False
    Tsurf:  float = 20.5 # Initial surface T
    Tbotm:  float = 20.0 # Initial bottom T
    Fsolar: float = 0
    Qnet:   float = 0
    R0:     float = 0
    abs_profile_type: str = 'I' # Jerlov absorption profile type (I, IA, IB, II)

    # --- Derived quantities ---
    f:          float = field(init=False)
    nt:         int   = field(init=False)
    cdrho:      float = field(init=False)
    sqrt_cdrho: float = field(init=False)
    sqrt_dt:    float = field(init=False)

    # tensor kwargs
    farr_kwargs: Dict[str, Any] = field(init=False, repr=False)
    carr_kwargs: Dict[str, Any] = field(init=False, repr=False)
    
    def finalize(self) -> None:
        
        # --- Types kwargs ---
        dev = torch.device(self.device)
        self.farr_kwargs = {"dtype": torch.float32,   "device": dev}
        self.carr_kwargs = {"dtype": torch.complex64, "device": dev}
        
        # --- Derived scalars ---
        self.f = 2 * 7.2921e-5 * math.sin(math.radians(self.lat))
        self.cdrho = self.Cd * self.rho_a / self.rho_o
        self.sqrt_cdrho = self.cdrho ** 0.5
        self.sqrt_dt = self.dt ** 0.5
        self.nt = int(self.tend / self.dt)
        
        # --- TKE closure ---
        self.tke_rmxl_min = 1.0e-6 / (self.tke_ediff * self.tke_min ** 0.5)
        self.tke_ri_cri = 2.0 / (2.0 + self.tke_ediss / self.tke_ediff)