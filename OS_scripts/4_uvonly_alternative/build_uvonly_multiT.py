#!/usr/bin/env python3
"""
UV-only build variant for species that have a genuine multi-temperature
measurement series available (as opposed to build_uvonly_flat.py, which
always writes flat/T-independent data). Use this when a species has several
literature UV cross-section files at different temperatures from the same
study -- gives real (if coarse) temperature dependence instead of a flat
value repeated across all 19 TMOL columns.

Method: each TMOL grid point (100-4000K) is assigned the cross-section
spectrum from whichever MEASURED temperature is nearest -- clamped, never
extrapolated beyond the min/max measured T. Wherever the T-series has no
wavelength coverage at all (most literature series only cover a narrow
band), falls back to a flat widest-coverage single file if one exists.

How to pick a T-series + fallback for a new species: list every local file
for the species (`ls /groups/astro/tbalduin/crosssecs/MPI/<SPECIES>_*.txt`),
count temperatures and datapoints per file. Prefer a same-study series with
>=2 distinct temperatures for the T-series; among candidates, more
temperature points and more datapoints per point both count in its favor.
For the fallback, pick whichever single file (any T) has the widest
wavelength span and/or most datapoints, to fill gaps outside the T-series'
narrow window. If a species has no real multi-T series (every file is
~room temperature from different studies), just use build_flat_only()
instead of build_species() -- see the N2H4 example below.

Same wavenumber-order bug history as build_uvonly_flat.py applies here --
convert AA->cm-1 BEFORE sorting, never after (1/x reverses order).

Below is the actual selection this project used for C3H6/C4H4/C6H6/N2H4,
kept as a worked example -- replace with your own species' file lists and
selections when reusing this script.
"""
import numpy as np
import astropy.units as u
from molmass import Formula
from scipy.interpolate import PchipInterpolator

CROSSSECS = '/groups/astro/tbalduin/crosssecs/MPI'
OUTDIR = '/groups/astro/tbalduin/marcs_opac_converter/combined_UV_IR/input_crossec_files/combined'
WN_GRID_FILE = '/groups/astro/tbalduin/python_scripts/data/MARCS_wavelengths/nwreal_7949.dat'

TMOL = [100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 500.0, 604.0,
        730.0, 882.0, 1065.0, 1286.0, 1554.0, 1878.0, 2269.0, 2741.0,
        3311.0, 4000.0]

wn_grid = 1e8 / np.loadtxt(WN_GRID_FILE)  # AA -> cm-1, same convention as the flat-builder script
wn_grid = np.sort(wn_grid)   # sort AFTER converting -- wavelength ascending inverts to wavenumber descending


def load_mpi(path):
    wl, xs = [], []
    with open(path, errors='replace') as f:
        for line in f:
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                w, x = float(parts[0]), float(parts[1])
            except ValueError:
                continue
            wl.append(w); xs.append(x)
    wl = np.asarray(wl); xs = np.asarray(xs)
    order = np.argsort(wl)
    return wl[order], xs[order]  # nm, cm^2/molecule


def spectrum_on_grid(path):
    """Return cross section (cm^2/molec) on wn_grid, 0 outside this file's own range."""
    wl_nm, xs = load_mpi(path)
    wn = 1e7 / wl_nm[::-1]
    xs = xs[::-1]
    interp = PchipInterpolator(wn, xs, extrapolate=False)
    out = interp(wn_grid)
    return np.where(np.isfinite(out), out, 0.0)


def write_os_file(molid, cs_per_T_g):
    """cs_per_T_g: (ktemp, nwnos) array in cm^2/g."""
    ktemp = len(TMOL)
    nwnos = len(wn_grid)
    out_path = f'{OUTDIR}/{molid}_UVonly_7949.dat'

    def fmt_tmol_lines(tmol):
        lines = []
        for i in range(0, len(tmol), 4):
            chunk = tmol[i:i + 4]
            lines.append(' '.join(f'{t:.1f} ,' for t in chunk))
        return lines

    reliso_extra = ' '.join(['0.0 ,'] * 14)
    with open(out_path, 'w') as fh:
        fh.write('&INPUTOSMOL \n')
        fh.write(f'MOLID = {molid} , \n')
        fh.write(f'KTEMP = {ktemp}, \n')
        tmol_lines = fmt_tmol_lines(TMOL)
        fh.write('TMOL = ' + tmol_lines[0] + ' \n')
        for tl in tmol_lines[1:]:
            fh.write(tl + ' \n')
        fh.write(f'NWNOS = {nwnos}, \n')
        fh.write('VKMS = 3.0 , \n')
        fh.write('KISO = 1, \n')
        fh.write(f'RELISO = 1.0 , {reliso_extra} \n')
        fh.write('L_PER_STELLAR = 0, \n')
        fh.write('LCHROM = 0 \n')
        fh.write('/ \n')
        for i, wn in enumerate(wn_grid):
            if wn == 0.0:
                continue
            cs_str = ' '.join(f'{cs_per_T_g[k, i]:.3e}' for k in range(ktemp))
            fh.write(f'{wn:.3f} {cs_str}\n')
    print(f'{molid}: wrote {out_path}')


def molec_to_g(molid, cs_molec):
    molmass = Formula(molid).mass * u.u
    fact = (u.cm**2 / molmass).to(u.cm**2 / u.g).value
    return cs_molec * fact


def build_species(molid, tseries_files_T, fallback_file):
    """tseries_files_T: list of (T_kelvin, filepath), sorted or not."""
    tseries_files_T = sorted(tseries_files_T, key=lambda x: x[0])
    Ts = np.array([t for t, _ in tseries_files_T])
    spectra_molec = [spectrum_on_grid(p) for _, p in tseries_files_T]  # each (nwnos,)

    fallback_molec = spectrum_on_grid(fallback_file) if fallback_file else np.zeros_like(wn_grid)

    ktemp = len(TMOL)
    nwnos = len(wn_grid)
    cs_g = np.zeros((ktemp, nwnos))

    fallback_g = molec_to_g(molid, fallback_molec)
    tseries_g = [molec_to_g(molid, s) for s in spectra_molec]

    for k, T in enumerate(TMOL):
        # clamp to nearest measured T (never extrapolate)
        idx = int(np.argmin(np.abs(Ts - T)))
        tcurve = tseries_g[idx]
        # where the T-series has no data (0), fall back to the flat coverage file
        col = np.where(tcurve > 0, tcurve, fallback_g)
        cs_g[k, :] = col

    write_os_file(molid, cs_g)
    n_tseries_pts = int((np.array(tseries_g) > 0).any(axis=0).sum())
    has_tseries = (np.array(tseries_g) > 0).any(axis=0)
    n_fallback_only = int((~has_tseries & (fallback_g > 0)).sum()) if fallback_file else 0
    print(f'  T-series measured T: {list(Ts)} K (clamped outside this range)')
    print(f'  grid points with real T-dependence: {n_tseries_pts}, '
          f'flat-fallback-only points: {n_fallback_only}, total nonzero: '
          f'{int((cs_g.max(axis=0) > 0).sum())} / {nwnos}')


def build_flat_only(molid, filepath):
    cs_molec = spectrum_on_grid(filepath)
    cs_g = molec_to_g(molid, cs_molec)
    ktemp = len(TMOL)
    full = np.tile(cs_g, (ktemp, 1))
    write_os_file(molid, full)
    print(f'  flat (T-independent), nonzero: {int((cs_g > 0).sum())} / {len(wn_grid)}')


# ============================================================================
# Worked example: the actual C3H6/C4H4/C6H6/N2H4 selection used in this
# project. Replace with your own species below.
# ============================================================================

# --- C3H6 ---
c3h6_tseries = [
    (223, f'{CROSSSECS}/C3H6_FahrNayak(1996)_223K_160-200nm.txt'),
    (233, f'{CROSSSECS}/C3H6_FahrNayak(1996)_233K_160-200nm.txt'),
    (253, f'{CROSSSECS}/C3H6_FahrNayak(1996)_253K_160-200nm.txt'),
    (273, f'{CROSSSECS}/C3H6_FahrNayak(1996)_273K_160-200nm.txt'),
    (295, f'{CROSSSECS}/C3H6_FahrNayak(1996)_295K_160-200nm.txt'),
    (313, f'{CROSSSECS}/C3H6_FahrNayak(1996)_313K_160-200nm.txt'),
    (333, f'{CROSSSECS}/C3H6_FahrNayak(1996)_333K_160-200nm.txt'),
]
print('=== C3H6 ===')
build_species('C3H6', c3h6_tseries, f'{CROSSSECS}/C3H6_Christianson(2021)_323K_125-240nm.txt')

# --- C4H4 ---
c4h4_tseries = [
    (223, f'{CROSSSECS}/C4H4_FahrNayak(1996)_223K_160-240nm.txt'),
    (233, f'{CROSSSECS}/C4H4_FahrNayak(1996)_233K_160-240nm.txt'),
    (253, f'{CROSSSECS}/C4H4_FahrNayak(1996)_253K_160-240nm.txt'),
    (273, f'{CROSSSECS}/C4H4_FahrNayak(1996)_273K_160-240nm.txt'),
    (295, f'{CROSSSECS}/C4H4_FahrNayak(1996)_295K_160-240nm.txt'),
    (313, f'{CROSSSECS}/C4H4_FahrNayak(1996)_313K_160-240nm.txt'),
    (333, f'{CROSSSECS}/C4H4_FahrNayak(1996)_333K_160-240nm.txt'),
]
print('=== C4H4 ===')
build_species('C4H4', c4h4_tseries, None)

# --- C6H6 ---
c6h6_tseries = [
    (253, f'{CROSSSECS}/C6H6_Fally(2009)_253K_239.236-270.269nm.txt'),
    (263, f'{CROSSSECS}/C6H6_Fally(2009)_263K_242.770-270.269nm.txt'),
    (273, f'{CROSSSECS}/C6H6_Fally(2009)_273K_239.236-270.269nm.txt'),
    (283, f'{CROSSSECS}/C6H6_Fally(2009)_283K_239.236-270.269nm.txt'),
    (293, f'{CROSSSECS}/C6H6_Fally(2009)_293K_239.236-270.269nm.txt'),
]
print('=== C6H6 ===')
build_species('C6H6', c6h6_tseries, f'{CROSSSECS}/C6H6_Dawes(2017)_298K_115-330nm.txt')

# --- N2H4 (no real multi-T series available; flat best-coverage file) ---
print('=== N2H4 ===')
build_flat_only('N2H4', f'{CROSSSECS}/N2H4_Vaghjiani(1993)_296K_191-291nm.txt')
