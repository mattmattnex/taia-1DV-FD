import ctypes as c

# Types for setting attributes and return types
ptr_int = c.POINTER(c.c_int)
ptr_dbl = c.POINTER(c.c_double)

# Each function from self.lib (libeco3m.so)
Eco3M_function_argtypes = {
    'eco3m_c_initialize':     [],
    'eco3m_c_step':           [],

    'eco3m_c_set_dims':       [c.c_int] * 5,
    'eco3m_c_set_bathymetry': [c.c_double],
    'eco3m_c_set_dt_bio':     [c.c_double],
    'eco3m_c_set_param':      [c.c_int, c.c_int, c.c_double, ptr_int],
    'eco3m_c_set_physics':    [ptr_dbl, ptr_dbl, c.c_double, c.c_double],
    'eco3m_c_set_dz':         [ptr_dbl],
    'eco3m_c_set_depth':      [ptr_dbl],
    'eco3m_c_set_conc':       [ptr_dbl],
    
    'eco3m_c_get_pos_nzmax':  [c.c_char_p, c.c_int],
    'eco3m_c_get_varname':    [c.c_int, c.c_char_p, c.c_char_p, c.c_char_p, c.c_int],
    'eco3m_c_get_nparam':     [c.c_int, ptr_int],
    'eco3m_c_get_flux':       [c.c_int, c.c_int, ptr_dbl],
    'eco3m_c_get_param':      [c.c_int, c.c_int, ptr_dbl, ptr_int],
    'eco3m_c_get_epar':       [ptr_dbl],
    'eco3m_c_get_dt_bio':     [ptr_dbl],
    'eco3m_c_get_nbvar':      [ptr_int],
    'eco3m_c_get_nzmax':      [ptr_int],
    'eco3m_c_get_nflux_par':  [ptr_int],
    'eco3m_c_get_conc':       [ptr_dbl],
    'eco3m_c_get_tend':       [ptr_dbl],
    'eco3m_c_get_eparz':      [ptr_dbl],
}