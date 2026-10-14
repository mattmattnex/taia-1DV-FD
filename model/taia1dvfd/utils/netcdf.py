import torch
from netCDF4 import Dataset
from pathlib import Path

from taia1dvfd.core import Model, State

class SaveNetCDF:
    """
    Ensemble statistics (mean, variance) written to NetCDF.

    Statistics are computed on the GPU and stored in a GPU buffer; the buffer is
    copied to the CPU and written to disk only every `flush_every` records.
    -> one GPU/CPU synchronisation and one NetCDF write per variable per flush,
       instead of several per variable per time step.

    Config (optional):
      save_every  : save one record every `save_every` calls to save()   (default 1)
      flush_every : number of records kept in the GPU buffer before writing (default 100)
    """
    def __init__(self, model: Model):
        self.cfg  = model.config
        self.grid = model.grid

        self.file = self.cfg.save_folder + self.cfg.file_name

        self.save_every  = self.cfg.save_every
        self.flush_every = self.cfg.flush_every

        self.k       = 0    # number of calls to save()
        self.n_saved = 0    # records already written to disk
        self.n_buf   = 0    # records currently in the GPU buffer

        # Will be filled in init_file:
        self.ds       = None
        self.var_time = None
        self.buf_time = None
        self.fields   = {}  # name -> metadata, nc variables and GPU buffer

    # =====================================================================
    #  File creation
    # =====================================================================
    def init_file(self, state: State):
        cfg = self.cfg
        fkw = cfg.farr_kwargs
        nz  = cfg.nz

        # --- Create a path ---
        path = Path(self.file)
        path.parent.mkdir(parents=True, exist_ok=True)

        # --- Open NetCDF file ---
        self.ds = Dataset(path, "w")

        # --- Dimensions ---
        self.ds.createDimension("time", None)
        self.ds.createDimension("zc", nz)          # cell centres
        self.ds.createDimension("zw", nz + 1)      # faces

        # --- Coordinates ---
        self.var_time = self.ds.createVariable("time", "f8", ("time",),
                                               chunksizes=(self.flush_every,))
        self.var_time.units = "seconds since start"

        var_zc = self.ds.createVariable("zc", "f4", ("zc",))
        var_zc[:] = self.grid.z_c.squeeze().cpu().numpy()
        var_zc.units = "m"

        var_zw = self.ds.createVariable("zw", "f4", ("zw",))
        var_zw[:] = self.grid.z_w.squeeze().cpu().numpy()
        var_zw.units = "m"

        # --- Time buffer (CPU, no GPU sync needed) ---
        self.buf_time = [0.0] * self.flush_every

        # --- Per-variable NetCDF variables and GPU buffers ---
        for name in cfg.var_to_save:
            tensor     = getattr(state, name)
            is_complex = torch.is_complex(tensor)

            # Vertical location from dimension 1
            z_location = None
            if tensor.dim() >= 2 and tensor.shape[1] == nz:
                z_location = "zc"
            elif tensor.dim() >= 2 and tensor.shape[1] == nz + 1:
                z_location = "zw"

            has_z = z_location is not None
            n_loc = tensor.shape[1] if has_z else 1

            if has_z:
                dims   = ("time", z_location)
                chunks = (self.flush_every, n_loc)
            else:
                dims   = ("time",)
                chunks = (self.flush_every,)

            # Statistics stored for this variable (order = order in the buffer)
            if is_complex:
                keys = ["mean_re", "mean_im", "var_re", "var_im"]
            else:
                keys = ["mean", "var"]

            nc_vars = {}
            for key in keys:
                nc_vars[key] = self.ds.createVariable(f"{name}_{key}", "f4", dims,
                                                      chunksizes=chunks)

            # GPU buffer: (flush_every, n_stats, n_loc)
            buffer = torch.zeros(self.flush_every, len(keys), n_loc, **fkw)

            self.fields[name] = {
                "is_complex": is_complex,
                "has_z":      has_z,
                "keys":       keys,
                "vars":       nc_vars,
                "buf":        buffer,
            }

        return state

    # =====================================================================
    #  Saving
    # =====================================================================
    def save(self, state: State):
        # Sub-sampling in time
        if self.k % self.save_every != 0:
            self.k += 1
            return

        # Time computed on the CPU (no sync with the GPU)
        self.buf_time[self.n_buf] = self.k * self.cfg.dt

        # Statistics computed and stored on the GPU (no sync)
        for name, meta in self.fields.items():
            tensor = getattr(state, name)
            meta["buf"][self.n_buf] = self._stats(tensor, meta["is_complex"])

        self.k     += 1
        self.n_buf += 1

        if self.n_buf == self.flush_every:
            self.flush()

    @staticmethod
    def _stats(tensor: torch.Tensor, is_complex: bool) -> torch.Tensor:
        """Ensemble statistics, shape (n_stats, n_loc). Stays on the GPU."""
        t = tensor.reshape(tensor.shape[0], -1)          # (ne, n_loc)

        if is_complex:
            mean = t.mean(dim=0)
            return torch.stack([mean.real,
                                mean.imag,
                                t.real.var(dim=0, correction=1 if t.shape[0] > 1 else 0),
                                t.imag.var(dim=0, correction=1 if t.shape[0] > 1 else 0)])
        else:
            return torch.stack([t.mean(dim=0),
                                t.var(dim=0, correction=1 if t.shape[0] > 1 else 0)])

    def flush(self):
        """Copy the buffered records to the CPU and write them in one block."""
        if self.n_buf == 0:
            return

        n  = self.n_buf
        sl = slice(self.n_saved, self.n_saved + n)

        # --- Time ---
        self.var_time[sl] = self.buf_time[:n]

        # --- Variables: one GPU -> CPU transfer per variable ---
        for name, meta in self.fields.items():
            block = meta["buf"][:n].cpu().numpy()        # (n, n_stats, n_loc)

            for i_stat, key in enumerate(meta["keys"]):
                if meta["has_z"]:
                    meta["vars"][key][sl, :] = block[:, i_stat, :]
                else:
                    meta["vars"][key][sl]    = block[:, i_stat, 0]

        self.n_saved += n
        self.n_buf    = 0
        self.ds.sync()                                   # data safe on disk if the run crashes later

    def close(self):
        if self.ds is not None:
            self.flush()
            self.ds.close()
            self.ds = None