import torch
from netCDF4 import Dataset
from pathlib import Path

from taia1dvfd.core import Model, State

class SaveNetCDF:
    def __init__(self, model: Model):
        self.cfg  = model.config
        self.grid = model.grid

        self.file = self.cfg.save_folder + self.cfg.file_name

        self.k = 0          # total step count
        self.n_saved = 0    # number of output records written

        # Will be filled in init_file:
        self.ds = None
        self.var_time = None
        self.var_zc = None
        self.var_zw = None
        self.fields = {}    # mapping name -> dict of nc variables

    def init_file(self, state: State):
        # --- Create a path ---
        path = Path(self.file)
        path.parent.mkdir(parents=True, exist_ok=True)

        # --- Open NetCDF file ---
        self.ds = Dataset(path, "w")

        # --- Dimensions ---
        self.ds.createDimension("time", None)
        self.ds.createDimension("zc", self.cfg.nz)          # center points
        self.ds.createDimension("zw", self.cfg.nz + 1)      # w-points (edges)

        # --- Coordinates ---
        self.var_time = self.ds.createVariable("time", "f8", ("time",))
        self.var_time.units = "seconds since start"

        self.var_zc = self.ds.createVariable("zc", "f4", ("zc",))
        self.var_zc[:] = self.grid.z_c.squeeze().detach().cpu().numpy()
        self.var_zc.units = "m"

        self.var_zw = self.ds.createVariable("zw", "f4", ("zw",))
        self.var_zw[:] = self.grid.z_w.squeeze().detach().cpu().numpy()
        self.var_zw.units = "m"

        # --- Prepare per-variable NetCDF vars for mean/var ---
        for name in self.cfg.var_to_save:
            tensor = getattr(state, name)
            is_complex = torch.is_complex(tensor)
            
            # Determine grid location (center or w-point)
            # Check dimension 1 of tensor
            has_z = tensor.dim() >= 2
            if has_z:
                z_dim = tensor.shape[1]
                if z_dim == self.cfg.nz:
                    z_location = "zc"
                elif z_dim == self.cfg.nz + 1:
                    z_location = "zw"
                else:
                    # Fallback: treat as no-z if dimension doesn't match
                    has_z = False
            
            entry = {}

            if has_z:
                dims = ("time", z_location)
            else:
                # Treat as scalar in z (e.g. ua)
                dims = ("time",)

            if is_complex:
                # mean(real), mean(imag), var(real), var(imag)
                entry["mean_re"] = self.ds.createVariable(f"{name}_mean_re", "f4", dims)
                entry["mean_im"] = self.ds.createVariable(f"{name}_mean_im", "f4", dims)
                entry["var_re"]  = self.ds.createVariable(f"{name}_var_re",  "f4", dims)
                entry["var_im"]  = self.ds.createVariable(f"{name}_var_im",  "f4", dims)
            else:
                entry["mean"] = self.ds.createVariable(f"{name}_mean", "f4", dims)
                entry["var"]  = self.ds.createVariable(f"{name}_var",  "f4", dims)

            self.fields[name] = {
                "is_complex": is_complex,
                "has_z": has_z,
                "z_location": z_location if has_z else None,
                "vars": entry,
            }

        return state

    def save(self, state: State):
        # Advance global step counter
        self.k += 1

        # Time index in NetCDF
        t_idx = self.n_saved
        self.var_time[t_idx] = float(state.time)

        # Loop over all variables we decided to save
        for name, meta in self.fields.items():
            tensor = getattr(state, name)
            is_complex = meta["is_complex"]
            has_z = meta["has_z"]
            ncvars = meta["vars"]

            t = tensor

            if has_z:
                # Expect shape (ne, nz, 1) or (ne, nz) or (ne, nz+1, 1) or (ne, nz+1)
                if t.dim() == 3 and t.shape[-1] == 1:
                    t = t.squeeze(-1)  # (ne, nz) or (ne, nz+1)

                if is_complex:
                    # Mean of full complex, var of real/imag separately
                    mean   = t.mean(dim=0)           # complex (nz,) or (nz+1,)
                    var_re = t.real.var(dim=0)       # real (nz,) or (nz+1,)
                    var_im = t.imag.var(dim=0)       # real (nz,) or (nz+1,)

                    mean_re   = mean.real.detach().cpu().numpy().astype("float32")
                    mean_im   = mean.imag.detach().cpu().numpy().astype("float32")
                    var_re_np = var_re.detach().cpu().numpy().astype("float32")
                    var_im_np = var_im.detach().cpu().numpy().astype("float32")

                    ncvars["mean_re"][t_idx, :] = mean_re
                    ncvars["mean_im"][t_idx, :] = mean_im
                    ncvars["var_re"][t_idx, :]  = var_re_np
                    ncvars["var_im"][t_idx, :]  = var_im_np

                else:
                    # Real-valued with z dimension
                    mean = t.mean(dim=0)
                    var  = t.var(dim=0)

                    mean_np = mean.detach().cpu().numpy().astype("float32")
                    var_np  = var.detach().cpu().numpy().astype("float32")

                    ncvars["mean"][t_idx, :] = mean_np
                    ncvars["var"][t_idx, :]  = var_np

            else:
                # No z-dimension: treat as scalar-ish (e.g. ua)
                if is_complex:
                    mean   = t.mean()              # complex scalar
                    var_re = t.real.var()          # real scalar
                    var_im = t.imag.var()          # real scalar

                    mean_re   = float(mean.real.detach().cpu().numpy())
                    mean_im   = float(mean.imag.detach().cpu().numpy())
                    var_re_np = float(var_re.detach().cpu().numpy())
                    var_im_np = float(var_im.detach().cpu().numpy())

                    ncvars["mean_re"][t_idx] = mean_re
                    ncvars["mean_im"][t_idx] = mean_im
                    ncvars["var_re"][t_idx]  = var_re_np
                    ncvars["var_im"][t_idx]  = var_im_np

                else:
                    mean = t.mean()
                    var  = t.var()

                    mean_np = float(mean.detach().cpu().numpy())
                    var_np  = float(var.detach().cpu().numpy())

                    ncvars["mean"][t_idx] = mean_np
                    ncvars["var"][t_idx]  = var_np

        # One more record written
        self.n_saved += 1

    def close(self):
        if self.ds is not None:
            self.ds.close()
            self.ds = None
