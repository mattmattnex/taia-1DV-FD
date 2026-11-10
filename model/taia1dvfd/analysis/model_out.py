import netCDF4 as nc

def get_model_variable(file_name:str, var_name: str, cplx: bool = False, kt:int | None = None):
        with nc.Dataset(file_name) as data:
            # --- Variables ---
            if cplx:
                var_m = data[f'{var_name}_mean_re'][:] + 1j*data[f'{var_name}_mean_im'][:]
                var_v = data[f'{var_name}_var_re'][:]  + 1j*data[f'{var_name}_var_im'][:]
            else:
                var_m = data[f'{var_name}_mean'][:]
                var_v = data[f'{var_name}_var'][:]

            # --- Coordinates ---
            z_c = data['zc'][:]
            z_w = data['zw'][:]
            t   = data['time'][:]

            # --- Find shape ---
            if len(var_m.shape) > 1:
                z = z_c if var_m.shape[1] == z_c.shape[0] else z_w
            else:
                z = 0.0
            
        if kt is not None:
            return z, var_m[kt, :], var_v[kt, :]
        else:
            return z, t, var_m, var_v