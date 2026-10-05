import ctypes
import os
import numpy as np

from .eco3m_cfunc import Eco3M_function_argtypes
from .eco3m_utils import _chdir, _f64, _ptr

_NBUF = 64  # scratch size for name strings

class Eco3MCoupler:
    """
    ctypes wrapper around libeco3m.so.

    Array convention
    ----------------
    Every 2D array crossing the boundary is (nz, nbvar) in FORTRAN order
    (first index varies fastest), matching TEND's native memory layout.

    Working directory
    -----------------
    eco3m_init_config resolves paths against eco3m_root_dir from the
    namelist (default "../") and writes OUTPUTS/<run_name>.log. Run from the
    MAKE/ directory, or pass workdir= to Eco3M().

    """

    def __init__(self, libpath="./libeco3m.so", workdir=None):
        self.workdir = workdir
        self.lib = ctypes.CDLL(os.path.abspath(libpath))
        self.nz = None
        self.nbvar = None
        self.varnames = []
        self.index = {}
        
        # Setting up argtypes for C bindings function (This is a header)
        for fn, args in Eco3M_function_argtypes.items():
            getattr(self.lib, fn).argtypes = args
            getattr(self.lib, fn).restype = None
            
        # Data container
        self.CONC = None

    def initialize(self, nz, nx_min=1, nx_max=1, ny_min=1, ny_max=1):
        """Initialize dimensions, read the namelist and config.ini, allocate everything."""
        
        # Initialize dimensions
        self.lib.eco3m_c_set_dims(nx_min, nx_max, ny_min, ny_max, nz)
        
        with _chdir(self.workdir):
            self.lib.eco3m_c_initialize()

        self.nz = self._int_out("eco3m_c_get_nzmax")
        self.nbvar = self._int_out("eco3m_c_get_nbvar")

        if self.nz <= 0 or self.nbvar <= 0:
            raise RuntimeError(
                f"implausible dimensions after init (nz={self.nz}, "
                f"nbvar={self.nbvar}) -- was set_dims() called first?"
            )

        self._read_varnames()
        
        # Initialize concentrations
        self.CONC = np.zeros((self.nz, self.nbvar), order="F")

    def step(self, conc, temp, sal, irrad, time_s):
        """Convenience: push everything, run one step, return TEND."""
        
        # Setup the step
        self.set_physics(temp, sal, irrad, time_s)
        self.set_conc(conc)
        
        # Compute the step
        with _chdir(self.workdir):
            self.lib.eco3m_c_step()
        return self.get_tend()

    # ------------------------------------------------------------------
    # Setters
    # ------------------------------------------------------------------
    def set_dz(self, dz):
        """Layer thicknesses, length nz. MUST be positive.

        compute_eparz_rgb uses dz without abs(), unlike the non-RGB path,
        so negative thicknesses make light grow exponentially with depth.
        """
        a = _f64(dz, self.nz, "dz")
        if np.any(a <= 0):
            raise ValueError("dz must be strictly positive (RGB path lacks abs())")
        self.lib.eco3m_c_set_dz(_ptr(a))

    def set_depth(self, depth):
        a = _f64(depth, self.nz, "depth")
        self.lib.eco3m_c_set_depth(_ptr(a))

    def set_bathymetry(self, value):
        self.lib.eco3m_c_set_bathymetry(float(value))

    def set_dt_bio(self, dt):
        """Set Eco3M's dt_bio. Process functions use it in positivity
        limiters, so it must match the step you actually take."""
        self.lib.eco3m_c_set_dt_bio(float(dt))

    def set_conc(self, conc):
        """Push concentrations, shape (nz, nbvar)."""
        a = self._as_state(conc, "conc")
        self.lib.eco3m_c_set_conc(_ptr(a))

    def set_physics(self, temp, sal, irrad, time_s):
        """Push the physical environment for this step.

        irrad must be NET (albedo-corrected) surface shortwave in W/m2:
        under COUPL, Eco3M does not apply albedo, assuming the physical
        model already did. The namelist 'albedo' entry is inactive.
        """
        t = _f64(temp, self.nz, "temp")
        s = _f64(sal, self.nz, "sal")
        self.lib.eco3m_c_set_physics(_ptr(t), _ptr(s),
                                     float(irrad), float(time_s))

    # ------------------------------------------------------------------
    # Getters
    # ------------------------------------------------------------------
    def get_eparz(self):
        """PAR profile, length nz. Cheapest check on dz sign and grid
        orientation: it must decay monotonically away from the surface."""
        out = np.zeros(self.nz, dtype=np.float64)
        self.lib.eco3m_c_get_eparz(_ptr(out))
        return out

    def get_epar(self):
        """Surface PAR just below the interface, after irr2par/irr_param."""
        return self._dbl_out("eco3m_c_get_epar")

    def get_flux(self, ivar, jvar):
        """One element of the flux matrix, 1-based, as a profile.

        FLUX_VAL(i,j) is the mass flux FROM variable i TO variable j.
        """
        out = np.zeros(self.nz, dtype=np.float64)
        self.lib.eco3m_c_get_flux(int(ivar), int(jvar), _ptr(out))
        return out
    
    def get_tend(self):
        out = np.zeros((self.nz, self.nbvar), dtype=np.float64, order="F")
        self.lib.eco3m_c_get_tend(_ptr(out))
        return out

    def get_conc(self):
        out = np.zeros((self.nz, self.nbvar), dtype=np.float64, order="F")
        self.lib.eco3m_c_get_conc(_ptr(out))
        return out

    # ------------------------------------------------------------------
    # parameters (calibration hook: no rebuild needed)
    # ------------------------------------------------------------------
    @property
    def dt_bio(self):
        """Biological time step from the namelist, in seconds."""
        return self._dbl_out("eco3m_c_get_dt_bio")

    @property
    def pos_nzmax(self):
        """'SURF' or 'BOTT': which end of the k index is the surface.

        'BOTT' means index nz_max is the bottom, so k=0 (Python) is the
        surface. Match your vertical grid to this.
        """
        buf = ctypes.create_string_buffer(_NBUF)
        self.lib.eco3m_c_get_pos_nzmax(buf, _NBUF)
        return buf.value.decode().strip()
    
    @property
    def n_flux_par(self):
        return self._int_out("eco3m_c_get_nflux_par")

    def n_param(self, iproc):
        return self._int_out("eco3m_c_get_nparam", int(iproc))

    def get_param(self, iproc, ipar):
        """Parameter value, both indices 1-based, matching call.inc."""
        v, ierr = ctypes.c_double(0.0), ctypes.c_int(0)
        self.lib.eco3m_c_get_param(int(iproc), int(ipar),
                                   ctypes.byref(v), ctypes.byref(ierr))
        if ierr.value != 0:
            raise IndexError(f"no parameter ({iproc}, {ipar})")
        return v.value

    def set_param(self, iproc, ipar, value):
        ierr = ctypes.c_int(0)
        self.lib.eco3m_c_set_param(int(iproc), int(ipar),
                                   float(value), ctypes.byref(ierr))
        if ierr.value != 0:
            raise IndexError(f"no parameter ({iproc}, {ipar})")

    def all_params(self):
        return {
            (ip, k): self.get_param(ip, k)
            for ip in range(1, self.n_flux_par + 1)
            for k in range(1, self.n_param(ip) + 1)
        }

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------
    def _int_out(self, fn, *args):
        v = ctypes.c_int(0)
        getattr(self.lib, fn)(*args, ctypes.byref(v))
        return v.value

    def _dbl_out(self, fn, *args):
        v = ctypes.c_double(0.0)
        getattr(self.lib, fn)(*args, ctypes.byref(v))
        return v.value
    
    def _read_varnames(self):
        """Query the name triplets and build the index map.

        Done at runtime rather than hardcoded: nbvar goes from 2 in the
        dummy config to ~40 in the full one, and a stale hardcoded index
        fails silently.
        """
        self.varnames, self.index = [], {}
        for ivar in range(1, self.nbvar + 1):
            bufs = [ctypes.create_string_buffer(_NBUF) for _ in range(3)]
            self.lib.eco3m_c_get_varname(ivar, *bufs, _NBUF)
            comp, scomp, elmt = (b.value.decode().strip() for b in bufs)
            name = f"{comp}_{scomp}_{elmt}"
            self.varnames.append(name)
            self.index[name] = ivar - 1  # 0-based for Python
        return self.varnames
    
    def _as_state(self, arr, name):
        """Coerce to an (nz, nbvar) float64 array in Fortran order."""
        a = np.asarray(arr, dtype=np.float64)
        if a.shape != (self.nz, self.nbvar):
            raise ValueError(
                f"{name}: expected shape ({self.nz}, {self.nbvar}), got {a.shape}"
            )
        return np.asfortranarray(a)