import torch
import netCDF4 as nc
from taia1dvfd.core import Model

def genererate_lotus_forcing(model: Model, window_len: int = 24):
    cfg = model.config
    
    # --- Parameters ---
    fkw = cfg.farr_kwargs
    ckw = cfg.carr_kwargs
    ne_model = cfg.ne
    dt_model = cfg.dt
    nt_model = cfg.nt
    wind_Ta = cfg.wind_Ta
    Ts_Ta   = cfg.wind_Ta

    # --- Lotus obs exctraction --- 
    LOTUS_MET_FILE = '/home/mnex/Desktop/1DV_FD/experiments/LOTUS/obs/lotusmet.epic'
    LOTUS_DT = 900
    LOTUS_NDAYS = 160
    LOTUS_NOBS = (LOTUS_NDAYS * 86400) // LOTUS_DT
    n_blocks = LOTUS_NOBS // window_len

    with nc.Dataset(LOTUS_MET_FILE, 'r') as f:
        f.set_auto_mask(False)
        LOTUS_UA = torch.tensor(f.variables['wnde'][:LOTUS_NOBS, 0, 0, 0], **fkw)
        LOTUS_VA = torch.tensor(f.variables['wndn'][:LOTUS_NOBS, 0, 0, 0], **fkw)
        LOTUS_TS = torch.tensor(f.variables['stmp'][:LOTUS_NOBS, 0, 0, 0], **fkw)
        
    bk_ua_mean, bk_ua_std, bk_ua_dt = _make_blocks(LOTUS_UA, LOTUS_DT, LOTUS_NOBS, window_len)
    bk_va_mean, bk_va_std, bk_va_dt = _make_blocks(LOTUS_VA, LOTUS_DT, LOTUS_NOBS, window_len)
    bk_Ts_mean, bk_Ts_std, bk_Ts_dt = _make_blocks(LOTUS_TS, LOTUS_DT, LOTUS_NOBS, window_len)
    
    t_forcing = torch.linspace(0, dt_model*nt_model, nt_model, **fkw) / 86400
    
    #--------------------------
    # --- Ornstein Ulenbeck ---
    #--------------------------
    
    # --- WIND ---
    c1_ua  = torch.exp(torch.tensor(-dt_model / wind_Ta))
    c2_ua  = torch.sqrt(1.0 - c1_ua**2)

    ua_forcing = torch.zeros(ne_model, nt_model, **ckw)
    ua_forcing[:, 0] = bk_ua_mean[0] + 1j*bk_va_mean[0]
    for kt in range(1, nt_model):
        # Index        
        t  = kt * dt_model
        kw = min(t // bk_ua_dt, n_blocks - 1)
        
        
        rd_cplx = torch.randn(ne_model, **ckw)
        ua_m   = bk_ua_mean[kw] + 1j*bk_va_mean[kw]
        ua_std = bk_ua_std[kw]  + 1j*bk_va_std[kw]
        ua_forcing[:, kt] = ua_m + c1_ua * (ua_forcing[:, kt-1] - ua_m) + c2_ua * ua_std * rd_cplx
        

    # --- SURFACE TEMP ---
    c1_Ts  = torch.exp(torch.tensor(-dt_model / Ts_Ta))
    c2_Ts  = torch.sqrt(1.0 - c1_Ts**2)
    
    Ts_forcing = torch.zeros(ne_model, nt_model, **fkw)
    Ts_forcing[:, 0] = bk_Ts_mean[0]
    for kt in range(1, nt_model):
        # Index        
        t  = kt * dt_model
        kw = min(t // bk_Ts_dt, n_blocks - 1)
        
        rd_float = torch.randn(ne_model, **fkw)
        Ts_m   = bk_Ts_mean[kw]
        Ts_std = bk_Ts_std[kw]
        
        Ts_forcing[:, kt] = Ts_m + c1_Ts * (Ts_forcing[:, kt-1] - Ts_m) + c2_Ts * Ts_std * rd_float

    return t_forcing, ua_forcing, Ts_forcing

def _make_blocks(var: torch.Tensor, dt_obs:float, n_obs:int, window_len:int):
    # ---Compute blocks ---
    n_blocks = n_obs // window_len
    blocks = var[:n_blocks * window_len].reshape(n_blocks, window_len)
    blocks_mean = blocks.mean(axis=1)
    blocks_std  = blocks.std(axis=1)
    blocks_dt = dt_obs * window_len
    
    return blocks_mean, blocks_std, blocks_dt