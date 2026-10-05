import os
import numpy as np
import ctypes as c
from contextlib import contextmanager

@contextmanager
def _chdir(path):
    """Temporarily change directory; Eco3M resolves file paths relatively."""
    if path is None:
        yield
        return
    prev = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev)

def _f64(arr, n, name):
    """Coerce to a contiguous float64 vector of length n."""
    a = np.ascontiguousarray(arr, dtype=np.float64).ravel()
    if a.size != n:
        raise ValueError(f"{name}: expected {n} values, got {a.size}")
    return a

def _ptr(a):
        return a.ctypes.data_as(c.POINTER(c.c_double))