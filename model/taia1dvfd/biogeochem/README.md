# `biogeochem` — Eco3M coupling

Drives the Fortran biogeochemical model **Eco3M** from the Python 1DV physics
model. Eco3M is compiled as a shared library and called once per timestep as a
pure function: given concentrations and a physical environment it returns
`TEND`, the source-minus-sink term **as a rate in per-second units**.

Transport, diffusion and time integration stay in Python. Eco3M's own TKE
physics is not compiled at all.

Current configuration: **Eco3M_MED_CN**, 15 state variables (copepods,
phytoplankton, heterotrophic bacteria, particulate and dissolved organic
matter, dissolved inorganic nitrogen).

---

## Files

| file | role |
|---|---|
| `bgc_coupling.py` | entry points used by the model: init and per-step |
| `eco3m_coupler.py` | `Eco3MCoupler`, the ctypes wrapper around `libeco3m.so` |
| `eco3m_cfunc.py` | ctypes signatures for every C entry point |
| `eco3m_utils.py` | `_chdir`, `_f64`, `_ptr` helpers |
| `eco3m_cini.py` | reads `C_INI` profiles, interpolates them onto the grid, validates quotas |
| `eco3m_test.py` | bring-up checks |

On the Fortran side:

| file | role |
|---|---|
| `Eco3M/SRC/src_python_coupl/eco3m_capi.F90` | the C-callable shim |
| `Eco3M/MAKE/Makefile.conf_PY`, `Makefile.calc_PY` | the shared-library build |
| `Eco3M/MAKE/libeco3m.so` | the built library |

---

## How the interface works

Eco3M is built with `-DECO3M_SUB`, which turns its main program into a
subroutine with no time loop. The whole body is:

```
eco3m_compute_eparz  ->  [compute_chl_c_ratio]  ->  eco3m_compute_sms
```

Python cannot call Eco3M directly: concentrations live in `VAR(:)%conc`, an
array of Fortran derived types whose memory layout is a compiler-internal
array descriptor. `eco3m_capi.F90` does the navigating and exposes flat
arrays through `bind(C)`, which also pins the exported symbol names so they
do not depend on the compiler.

Per step:

```
set_physics(temp, sal, irrad, time)   ->  TEMP_BIO, SAL_BIO, irrad, tps
set_conc(conc)                        ->  VAR(ivar)%conc, looped over ivar
eco3m_c_step()                        ->  eco3m_main
get_tend()                            ->  TEND, (nz, nbvar), per second
```

**Array convention:** every 2D array crossing the boundary is `(nz, nbvar)` in
**Fortran order**. Same shape in C order is a different arrangement of the same
bytes and gives scrambled profiles with no error, so the wrapper enforces it.

**Orientation:** `pos_nzmax = "BOTT"`, so index 0 is the surface on both sides.

---

## Running

```python
from taia1dvfd.biogeochem import bgc_coupling_init, bgc_coupling_step

bgc = bgc_coupling_init(model)          # dims, init, grid, initial conditions
...
bgc_coupling_step(model, bgc)           # once per timestep
```

The initialisation order is not negotiable:

```python
bgc.initialize(nz=cfg.nz)               # set_dims THEN eco3m_c_initialize
bgc.set_dt_bio(cfg.dt)                  # process functions read dt_bio
bgc.set_dz(np.full(cfg.nz, dz))         # strictly positive
bgc.CONC = load_cini(CINI, DEPTH, z, bgc.varnames)
```

Under `ECO3M_SUB` the namelist group carrying `nx_min..nz_max` is not read, and
`eco3m_init_config` allocates six arrays from those values. Calling
`eco3m_c_initialize` before `eco3m_c_set_dims` allocates on uninitialised
integers.

Eco3M resolves its own file paths relatively, so the wrapper `chdir`s into
`workdir` (which must be `Eco3M/MAKE`) around init and each step.

---

## Configuration

Three layers, all under `Eco3M/CONFIG_ECO3M/`:

```
PREPAR_CONFIG_BIO/config.ini_txt_*   flux matrix, symbolic parameter names
PREPAR_CONFIG_BIO/param.txt_*        parameter name -> value table
        | text2numb
        v
config.ini                           SAME matrix, numeric literals  (GENERATED)
modele.def                           process name -> Fortran function + arg list
F_PROCESS/**/*.inc                   the process function library
        | eco3m_ini.exe
        v
MAKE/call.inc                        the generated flux equations  (GENERATED)
```

Never hand-edit `config.ini` or `call.inc`.

### Changing a parameter value

Two options.

**At runtime, no rebuild** — parameters live in `FLUX_PAR` and are writable:

```python
bgc.get_param(iproc, ipar)              # both 1-based, matching call.inc
bgc.set_param(iproc, ipar, value)
bgc.all_params()
```

This is the calibration hook. Index order follows the order sub-processes
appear in the flux matrix; read `MAKE/call.inc` to map index to function.

**Permanently** — edit `param.txt_*`, then re-run stages A, B and C below.

### Changing the model structure

Variables, processes or the flux matrix: edit `config.ini_txt_*`, and if you
add a process function, put its `.inc` under `F_PROCESS/<CATEGORY>/`, register
it in `modele.def` as

```
shortname:fortran_function_name:nbpar:arg1:arg2:...
```

and **add an `#include` line to `SRC/src_eco3m/eco3m_fprocess.F90`**. That file
is a hand-maintained list — nothing generates it. A function registered in
`modele.def` but not included there fails at link time.

Note `modele.def` declares a maximum argument count; calling a function with
fewer is normal and the functions branch on what they receive. The
"Number of parameters ... is wrong" warnings at startup are expected.

### Initial conditions

`eco3m_cini.py` reads the `C_INI_DYFAMED` `.DAT` profiles and interpolates them
onto the model grid:

```python
conc = load_cini(cini_dir, depth_file, z_target, bgc.varnames)
```

`z_target` is **positive depths, surface first** — pass `-grid.z_c`, not
`grid.z_c`.

This is not a plain interpolation. Eco3M_MED_CN is a Droop quota model, so
cell abundance is interpolated directly while intracellular contents go
through their quotas (`C/cell`, `N/cell`, `Chl/C`) and are rebuilt as
`cell x quota`. Interpolating contents directly lets quotas drift outside
`[Qmin, Qmax]`, where the Droop expressions are meaningless.

The loader also works around two problems in the supplied `C_INI` set:

- `zoo_COP_N`, `phy_PHY_C`, `phy_PHY_N` are stale and are **regenerated** as
  `cell x constant quota`, with the fractions from `make_cini_1DV_CN_2026.f90`
  (carbon near `Qmin`, nitrogen near `Qmax` — a late-winter state).
  `phy_PHY_Chl` follows carbon.
- `mop_MOP_C/N` and `mod_MODLC_C`/`mod_MODLN_N` were overwritten with copies of
  bacterial biomass, so they are read from the `..._V0` library instead.

Pass `regenerate={}` or `overrides={}` to disable either.

Every build ends in `validate()`, which checks positive cell abundance, quotas
inside their bounds, C:N inside `[3, 22]` and Chl:C under `CHL_C_MAX`. It
raises rather than letting a bad state become a `NaN` later.

The `.DAT` files carry no depth coordinate — depths come from
`OUTPUTS_baseline_tke/sim_eco3m_1DV_depth.dat` (41 water levels, 0.5 m to
2228 m; level 42 is a NEMO seafloor cell and is dropped).

---

## Rebuilding Eco3M

```bash
cd Eco3M/CONFIG_ECO3M/PREPAR_CONFIG_BIO        # A: symbolic -> numeric
gfortran text2numb.f90 -o text2numb
./text2numb                                     # asks for the two filenames
cp config.ini_txt_conv ../config.ini

cd ../../MAKE                                   # B: generate call.inc
make -f Makefile.ini clean && make -f Makefile.ini
./eco3m_ini.exe

make -f Makefile.calc_PY clean && \
make -f Makefile.calc_PY MODE=debug             # C: build libeco3m.so
```

`MODE=debug` adds bounds checking and `-finit-integer`, which turns the
dimensions-unset trap into a loud allocation failure. Drop it for production.

| what changed | stages |
|---|---|
| parameter value | A, B, C (or `set_param` at runtime) |
| flux matrix, variables, processes | A, B, C |
| a `.inc` function body | C |
| Eco3M source or CPP keys | C |
| namelist | none — just re-copy it |

**The namelist lives in `MAKE/`.** Under `ECO3M_SUB` Eco3M looks for
`namelist_eco3m` in the working directory, not `../CONFIG_ECO3M/`. Two copies
exist; only the `MAKE/` one is read.

A stale `call.inc` against a new `config.ini` compiles cleanly and computes the
old biogeochemistry, so re-run stage B after any structural change.

---

## Gotchas

- **`dz` must be positive.** `compute_eparz_rgb` uses `dz` without `abs()`,
  unlike the non-RGB path; negative thicknesses make light grow with depth.
- **Irradiance must be albedo-corrected.** Under `COUPL`, Eco3M does not apply
  albedo — it assumes the physics already did. The namelist `albedo` is dead.
- **`fichirrad` must end in `CODEPHYS`**, or Eco3M computes its own light and
  silently ignores what Python passes.
- **`dt_bio` is not cosmetic.** The process functions read it in positivity
  limiters, so call `set_dt_bio()` with the step you actually take.
- **Quotas overshoot `Qmax` slightly** under forward Euler — about 0.2% at
  `dt = 240 s`, 0.04% at `dt = 60 s`, first-order in `dt`. Sub-step the biology
  if it matters; clamping breaks nitrogen conservation.
- **Sinking is not handled.** `mop_MOP_C/N` are particulate and should sink;
  Eco3M's `BIOSEDIM` lived in the discarded TKE code. Vertical export is the
  Python model's job.
- **`kextinc_rgb` is never set at `k = nz_max`** (`compute_rgb_kextinc` loops
  to `nz_max-1`). Harmless under `BOTT`, where that is the bottom cell.

### Conservation, as a test

No nitrogen variable has a diagonal entry in the flux matrix, so **nitrogen is
conserved exactly** across `zoo_COP_N`, `phy_PHY_N`, `bac_BAC_N`, `mop_MOP_N`,
`mod_MODLN_N`, `mid_DIN_N`. Carbon is not, because respiration and
photosynthesis exchange with CO2. Both are worth asserting — the first catches
index and ordering bugs, the second confirms the diagonal is firing.

### Open questions for the config author

- `mod_MODLC_C` and `mod_MODLN_N`: separate sub-compartments, but V0 has 12x
  more N than C. Intended?
- Initial `Chl_tot` is ~5e-9, below the floor of the `kRGB61.txt` lookup
  (`~1.06e-8`), so the chlorophyll-light feedback is currently clamped to the
  clearest-water row.