from __future__ import annotations

import os
import re

import numpy as np

# ----------------------------------------------------------------------
# Model structure. Eco3M_MED_CN, 15 state variables.
# ----------------------------------------------------------------------

# Quota keys match quotas_sim.txt: phytoplankton is PHY, zooplankton MZOO.
ORGANISMS = {
    "zoo_COP": {"quota_tag": "MZOO", "elements": ("C", "N")},
    "phy_PHY": {"quota_tag": "PHY", "elements": ("C", "N"), "chl": True},
    "bac_BAC": {"quota_tag": "BAC", "elements": ("C", "N")},
}

# No quota constraint: interpolate straight. f_hydrol's Qmin/Qmax arguments
# are the BACTERIAL quotas, not the detritus's, so mop/mod C:N is free.
DIRECT_VARS = (
    "mop_MOP_C",
    "mop_MOP_N",
    "mod_MODLC_C",
    "mod_MODLN_N",
    "mid_DIN_N",
)

# Quota = w_min*Qmin + w_max*Qmax, from make_cini_1DV_CN_2026.f90 with the
# PHYL -> PHY names corrected. Weight on Qmax sets the physiological state:
# 0.10 is near-limiting, 0.75 near-replete.
#     varname -> (quota_tag, element, w_min, w_max)
REGEN_FRACTIONS = {
    "zoo_COP_N": ("MZOO", "N", 0.25, 0.75),
    "phy_PHY_C": ("PHY", "C", 0.90, 0.10),
    "phy_PHY_N": ("PHY", "N", 0.25, 0.75),
}

# ChlCmoy from the generator: C:Chl = 70 g C / g Chl, times 12 g C per mol C.
# Units g Chl / mol C, matching config.ini. Applied when phy carbon is
# regenerated, so chlorophyll stays consistent with it.
CHL_C_REGEN = 12.0 / 70.0

# The pristine library the Readme calls BIBLIOTHEQUE_CINI/CINI_V0, nested
# inside C_INI_DYFAMED.
V0_SUBDIR = "C_INI_DYFAMED_ECO3M_MED_EMCCC_V0"

# make_cini lines 85 and 91 never reassign ficell for gf==4 and gf==5, so
# nom stays "bac" from gf==3. Those iterations take the living-organism
# branch and write detritus as (bacterial cell) x (bacterial quota). The
# four files below therefore hold identical copies of bacterial biomass;
# the V0 originals are used instead.
GENERATOR_CLOBBERED = ("mop_MOP_C", "mop_MOP_N", "mod_MODLC_C", "mod_MODLN_N")

# NEMO seafloor cell -- not ocean.
SEAFLOOR_LEVEL = 42

# From the f_rchl(...) call in config.ini, argument 2. Units: g Chl / mol C.
CHL_C_MAX_DEFAULT = 1.20


# ----------------------------------------------------------------------
# File readers
# ----------------------------------------------------------------------

def read_dat(path: str, drop_seafloor: bool = True) -> np.ndarray:
    """Read a C_INI .DAT profile ("i j k value") into a 1D array by level.

    Returns values ordered by k ascending, i.e. surface first.
    """
    levels, values = [], []
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                k = int(parts[2])
                v = float(parts[3].replace("D", "E").replace("d", "e"))
            except ValueError:
                continue
            if drop_seafloor and k == SEAFLOOR_LEVEL:
                continue
            levels.append(k)
            values.append(v)

    if not values:
        raise ValueError(f"{path}: no data rows parsed")

    order = np.argsort(levels)
    return np.asarray(values, dtype=np.float64)[order]


def read_depths(path: str, drop_seafloor: bool = True) -> np.ndarray:
    """Read sim_eco3m_1DV_depth.dat into positive depths, surface first."""
    levels, depths = [], []
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                k = int(parts[0])
                d = float(parts[1])
            except ValueError:
                continue            # header line
            if drop_seafloor and k == SEAFLOOR_LEVEL:
                continue
            levels.append(k)
            depths.append(d)

    order = np.argsort(levels)
    z = np.asarray(depths, dtype=np.float64)[order]

    if not np.all(np.diff(z) > 0):
        raise ValueError(f"{path}: depths are not strictly increasing")
    return z


def read_quotas(path: str) -> dict[str, float]:
    """Parse quotas_sim.txt (KEY=VALUE lines; a Fortran include snippet)."""
    out = {}
    pat = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*([0-9eEdD.+-]+)")
    with open(path) as fh:
        for line in fh:
            m = pat.match(line)
            if m:
                key, raw = m.groups()
                out[key] = float(raw.replace("D", "E").replace("d", "e"))
    return out


def _dat_path(cini_dir: str, name: str, overrides: dict[str, str]) -> str:
    """Path to a variable's .DAT, honouring any per-variable directory."""
    return os.path.join(overrides.get(name, cini_dir), f"{name}.DAT")


def default_overrides(cini_dir: str) -> dict[str, str]:
    """Route the generator-clobbered variables to the V0 library, if present."""
    v0 = os.path.join(cini_dir, V0_SUBDIR)
    if not os.path.isdir(v0):
        return {}
    return {
        name: v0
        for name in GENERATOR_CLOBBERED
        if os.path.exists(os.path.join(v0, f"{name}.DAT"))
    }


def resolve_regen_quotas(
    quotas: dict[str, float],
    fractions: dict[str, tuple[str, str, float, float]],
) -> dict[str, float]:
    """Turn (tag, element, w_min, w_max) specs into constant quota values."""
    out = {}
    for name, (tag, elt, w_min, w_max) in fractions.items():
        kmin, kmax = f"Q_MIN_{tag}_{elt}", f"Q_MAX_{tag}_{elt}"
        if kmin not in quotas or kmax not in quotas:
            raise KeyError(
                f"{name}: quotas_sim.txt is missing {kmin} and/or {kmax}"
            )
        out[name] = w_min * quotas[kmin] + w_max * quotas[kmax]
    return out


# ----------------------------------------------------------------------
# Interpolation
# ----------------------------------------------------------------------

def _interp(z_src: np.ndarray, v_src: np.ndarray, z_dst: np.ndarray) -> np.ndarray:
    """Linear interpolation, clamped to the end values outside the source range.

    np.interp extends with the first/last value rather than extrapolating.
    That is what we want at the surface: at nz=512 the top cell centre
    (0.293 m) sits above the shallowest source level (0.506 m), and holding
    the surface value constant is right for a mixed layer.
    """
    return np.interp(z_dst, z_src, v_src)


def load_cini(
    cini_dir: str,
    depth_file: str,
    z_target: np.ndarray,
    varnames: list[str],
    chl_c_max: float = CHL_C_MAX_DEFAULT,
    quotas_file: str | None = None,
    regenerate: dict | None = None,
    overrides: dict[str, str] | None = None,
    verbose: bool = True,
) -> np.ndarray:
    """Build the (nz, nbvar) initial-condition array for Eco3M.

    Parameters
    ----------
    cini_dir   : directory holding the .DAT files
    depth_file : sim_eco3m_1DV_depth.dat
    z_target   : POSITIVE depths in metres, surface first (e.g. -grid.z_c)
    varnames   : bgc.varnames, which fixes the column order
    quotas_file: quotas_sim.txt; defaults to cini_dir/quotas_sim.txt
    regenerate : which contents to rebuild as cell x constant quota.
                 None  -> REGEN_FRACTIONS (the three stale profiles)
                 {}    -> none; use the source ratios as they are
                 dict  -> {varname: (tag, element, w_min, w_max)}
    overrides  : per-variable directory, for files to read elsewhere.
                 None -> the V0 library for GENERATOR_CLOBBERED
                 {}   -> read everything from cini_dir

    Returns
    -------
    (nz, nbvar) float64, Fortran-ordered, ready for bgc.set_conc().
    """
    z_target = np.ascontiguousarray(z_target, dtype=np.float64).ravel()
    if np.any(np.diff(z_target) <= 0):
        raise ValueError(
            "z_target must be strictly increasing (surface first, positive "
            "metres). For the torch grid, pass -grid.z_c, not grid.z_c."
        )

    if quotas_file is None:
        cand = os.path.join(cini_dir, "quotas_sim.txt")
        quotas_file = cand if os.path.exists(cand) else None

    quotas = read_quotas(quotas_file) if quotas_file else {}

    fractions = REGEN_FRACTIONS if regenerate is None else regenerate
    if fractions and not quotas:
        raise ValueError(
            "regeneration needs quotas_sim.txt; pass quotas_file= or "
            "regenerate={} to use the source ratios"
        )
    regen = resolve_regen_quotas(quotas, fractions) if fractions else {}

    if overrides is None:
        overrides = default_overrides(cini_dir)

    z_src = read_depths(depth_file)
    nz = z_target.size
    profiles: dict[str, np.ndarray] = {}

    if verbose:
        print(f"source grid : {z_src.size} levels, {z_src[0]:.3f} m to "
              f"{z_src[-1]:.1f} m (level {SEAFLOOR_LEVEL} dropped)")
        print(f"target grid : {nz} levels, {z_target[0]:.3f} m to "
              f"{z_target[-1]:.1f} m")
        if z_target[0] < z_src[0]:
            print(f"  note: top cell ({z_target[0]:.3f} m) is above the "
                  f"shallowest source level ({z_src[0]:.3f} m) "
                  f"-> clamped to the surface value")
        for name, q in sorted(regen.items()):
            print(f"  regenerating {name:<14} quota = {q:.4e} (constant)")
        for name, d in sorted(overrides.items()):
            print(f"  reading      {name:<14} from {os.path.basename(d)}/")

    # --- organisms: cell from file, contents via quota -------------------
    for org, spec in ORGANISMS.items():
        cell_name = f"{org}_cell"
        cell_src = read_dat(_dat_path(cini_dir, cell_name, overrides))

        if np.any(cell_src <= 0):
            raise ValueError(
                f"{cell_name}.DAT has non-positive values; the quota "
                f"functions divide by cell abundance"
            )

        cell_dst = _interp(z_src, cell_src, z_target)
        profiles[cell_name] = cell_dst

        contents = {}
        for elt in spec["elements"]:
            name = f"{org}_{elt}"
            if name in regen:
                quota_dst = np.full(nz, regen[name])          # constant
            else:
                src = read_dat(_dat_path(cini_dir, name, overrides))
                quota_dst = _interp(z_src, src / cell_src, z_target)
            contents[elt] = cell_dst * quota_dst
            profiles[name] = contents[elt]

        if spec.get("chl"):
            name = f"{org}_Chl"
            if f"{org}_C" in regen:
                # Carbon was rebuilt, so chlorophyll must follow it or the
                # Chl:C ratio inherits the stale carbon.
                ratio_dst = np.full(nz, CHL_C_REGEN)
                if verbose:
                    print(f"  regenerating {name:<14} Chl:C = "
                          f"{CHL_C_REGEN:.4f} g Chl/mol C (follows "
                          f"{org}_C)")
            else:
                src = read_dat(_dat_path(cini_dir, name, overrides))
                c_src = read_dat(_dat_path(cini_dir, f"{org}_C", overrides))
                ratio_dst = _interp(z_src, src / c_src, z_target)
            profiles[name] = contents["C"] * ratio_dst

    # --- everything else: straight interpolation --------------------------
    for name in DIRECT_VARS:
        src = read_dat(_dat_path(cini_dir, name, overrides))
        profiles[name] = _interp(z_src, src, z_target)

    # --- assemble in the model's own variable order -----------------------
    missing = [v for v in varnames if v not in profiles]
    if missing:
        raise KeyError(
            f"no profile built for {missing}. Known: {sorted(profiles)}. "
            f"If the model expects mid_NO3_N, this config file set is the "
            f"wrong one -- Eco3M_MED_CN uses mid_DIN_N."
        )

    conc = np.zeros((nz, len(varnames)), dtype=np.float64, order="F")
    for j, name in enumerate(varnames):
        conc[:, j] = profiles[name]

    if quotas:
        validate(conc, varnames, quotas, chl_c_max=chl_c_max, verbose=verbose)
    elif verbose:
        print("  warning: no quotas file, skipping quota checks")

    return conc


# ----------------------------------------------------------------------
# Pre-flight validation
# ----------------------------------------------------------------------

def validate(
    conc: np.ndarray,
    varnames: list[str],
    quotas: dict[str, float],
    chl_c_max: float = CHL_C_MAX_DEFAULT,
    verbose: bool = True,
) -> None:
    """Check the state is admissible before the first step().

    One assertion here beats hunting a NaN through 15 coupled variables.
    """
    idx = {n: i for i, n in enumerate(varnames)}
    problems: list[str] = []

    if not np.all(np.isfinite(conc)):
        problems.append("non-finite values in the initial state")
    if np.any(conc < 0):
        bad = [varnames[j] for j in range(conc.shape[1]) if conc[:, j].min() < 0]
        problems.append(f"negative concentrations in {bad}")

    for org, spec in ORGANISMS.items():
        tag = spec["quota_tag"]
        cell_name = f"{org}_cell"
        if cell_name not in idx:
            continue

        cell = conc[:, idx[cell_name]]
        if np.any(cell <= 0):
            problems.append(f"{cell_name}: non-positive cell abundance "
                            f"(min {cell.min():.3e}) -- quotas divide by this")
            continue

        for elt in spec["elements"]:
            name = f"{org}_{elt}"
            if name not in idx:
                continue
            qmin = quotas.get(f"Q_MIN_{tag}_{elt}")
            qmax = quotas.get(f"Q_MAX_{tag}_{elt}")
            if qmin is None or qmax is None:
                continue

            q = conc[:, idx[name]] / cell
            if q.min() < qmin or q.max() > qmax:
                problems.append(
                    f"{name}: quota [{q.min():.4e}, {q.max():.4e}] outside "
                    f"[{qmin:.4e}, {qmax:.4e}]"
                )
            elif verbose:
                frac = (q.max() - qmin) / (qmax - qmin)
                print(f"  {name:<16} quota {q.min():.4e}..{q.max():.4e}  "
                      f"({frac:6.1%} of range)")

        # molar C:N ratio
        cn_min = quotas.get(f"Q_MIN_{tag}_CN")
        cn_max = quotas.get(f"Q_MAX_{tag}_CN")
        cname, nname = f"{org}_C", f"{org}_N"
        if cn_min is not None and cname in idx and nname in idx:
            n = conc[:, idx[nname]]
            if np.any(n <= 0):
                problems.append(f"{nname}: non-positive, C:N undefined")
            else:
                cn = conc[:, idx[cname]] / n
                if cn.min() < cn_min or cn.max() > cn_max:
                    problems.append(
                        f"{org}: C:N [{cn.min():.3f}, {cn.max():.3f}] outside "
                        f"[{cn_min:.3f}, {cn_max:.3f}]"
                    )
                elif verbose:
                    print(f"  {org+' C:N':<16} {cn.min():.3f}..{cn.max():.3f}")

        if spec.get("chl"):
            chl_name = f"{org}_Chl"
            if chl_name in idx and cname in idx:
                r = conc[:, idx[chl_name]] / conc[:, idx[cname]]
                if r.max() > chl_c_max:
                    problems.append(
                        f"{chl_name}: Chl:C max {r.max():.4f} exceeds "
                        f"CHL_C_MAX {chl_c_max:.4f} g Chl/mol C"
                    )
                elif verbose:
                    print(f"  {org+' Chl:C':<16} {r.min():.4f}..{r.max():.4f} "
                          f"(max allowed {chl_c_max:.2f})")

    if problems:
        raise ValueError(
            "initial conditions are not admissible:\n  - "
            + "\n  - ".join(problems)
        )
    if verbose:
        print("  all quota and ratio checks passed")