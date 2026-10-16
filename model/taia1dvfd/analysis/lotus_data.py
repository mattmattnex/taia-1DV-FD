from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class LotusData:
    # ------------ Observations ------------
    # Velocity from Price et al. 1987
    CWC_mean = np.array([0.046,  0.028,  0.020,  0.004]) #, np.nan,  0.006,  0.011])    
    CWC_std  = np.array([0.012,  0.007,  0.007,  0.004]) #, np.nan,  0.003,  0.006])
    DWC_mean = np.array([0.010, -0.003, -0.002, -0.005]) #, np.nan, -0.002,  0.007])
    DWC_std  = np.array([0.007,  0.004,  0.005,  0.004]) #, np.nan,  0.004,  0.006])
    z         = np.array([-5.0,   -10.0,  -15.0,  -25.0]) #,  -50.0,  -75.0, -100.0])
    # -------------------------------------

    u =  DWC_mean[np.newaxis, :]
    v = -CWC_mean[np.newaxis, :]
    n =  z.shape[0]