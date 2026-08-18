#!/usr/bin/env python3
"""
Build UV-only MARCS OS intermediate files for species that have a UV cross
section on hand (Leiden/SWRI/MPI-Mainz format) but no rovibrational
linelist to run through os26 -- see ../README.md, "Path B" for when to use
this instead of the os26 pipeline in 2_os26_crosssections/.

Writes flat/T-independent files (UV absorption has no measured temperature
dependence in these sources) on the standard MARCS 7949-point wavenumber
grid, replicated across all 19 TMOL columns (100-4000K). For a species with
a genuine multi-temperature literature series available, use
build_uvonly_multiT.py instead to get real (clamped) temperature
dependence.

Edit SPECIES below and re-run for new species. `kind` selects the reader:
  leiden -- 4-column (or 3, if no photoionization channel) Leiden-format
            file: wavelength(nm), photoabsorption, photodissociation,
            [photoionization]. Cross section = max over all channel columns.
  swri   -- SWRI-format file: multiple reference/branching-ratio tables
            concatenated in one file, each preceded by free-text headers of
            varying length. Reads ONLY the first contiguous numeric block
            (the "Total" or sole first channel) and stops there -- do not
            let it run into a second table further down the file. Column 0
            is Lambda in Angstrom, column 1 is the total cross section
            (verified true for HNO4/N2O5/NO3/HNCO's SWRI files -- the
            DMS-era script's old default of column 3 was WRONG, don't reuse
            it blindly for a new species without checking that file's own
            column layout).
  mpi    -- plain 2-column MPI-Mainz format: wavelength(nm) <tab>
            cross-section(cm^2/molec), no header. Negative near-zero noise
            values at measurement edges get clipped to 0 same as
            everywhere else in this script.

Output: input_crossec_files/combined/<MOLID>_UVonly_7949.dat (namelist-
format OS text file) -- feed into 3_convert_opac/ next.

BUG HISTORY, read before touching the wn_grid conversion below: the
original version of this script did
    wn_grid = np.sort(np.loadtxt(WN_GRID_FILE)); wn_grid = 1e8/wn_grid
which sorts wavelength ascending THEN inverts to wavenumber -- since 1/x
reverses order, this silently wrote every one of 29 species' output files
with wavenumbers in descending order (MARCS requires ascending, refuses to
load anything else). Convert to wavenumber FIRST, sort second -- exactly as
done below. See UVonly_wavenumber_order_bug_report.md in this repo's
combined_UV_IR checkout for the full incident writeup if it exists, or just
don't reorder this again.
"""
import numpy as np
import astropy.units as u
from molmass import Formula
from scipy.interpolate import PchipInterpolator

CROSSSECS = '/groups/astro/tbalduin/crosssecs'
OUTDIR = '/groups/astro/tbalduin/marcs_opac_converter/combined_UV_IR/input_crossec_files/combined'
WN_GRID_FILE = '/groups/astro/tbalduin/python_scripts/data/MARCS_wavelengths/nwreal_7949.dat'

TMOL = [100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 500.0, 604.0,
        730.0, 882.0, 1065.0, 1286.0, 1554.0, 1878.0, 2269.0, 2741.0,
        3311.0, 4000.0]

SPECIES = [
    # (molid, kind, relpath under CROSSSECS)
    ('C2H',    'leiden', 'Leiden/C2H.txt'),
    ('C2H3',   'leiden', 'Leiden/C2H3.txt'),
    ('C2H5',   'leiden', 'Leiden/C2H5.txt'),
    ('C3H3',   'leiden', 'Leiden/C3H3.txt'),
    ('CH2',    'leiden', 'Leiden/CH2.txt'),
    ('CH3CHO', 'leiden', 'Leiden/CH3CHO.txt'),
    ('CH3NH2', 'leiden', 'Leiden/CH3NH2.txt'),
    ('CH3SH',  'leiden', 'Leiden/CH3SH.txt'),
    ('HCO',    'leiden', 'Leiden/HCO.txt'),
    ('HNC',    'leiden', 'Leiden/HNC.txt'),
    ('HNCO',   'leiden', 'Leiden/HNCO.txt'),
    ('NH2',    'leiden', 'Leiden/NH2.txt'),
    ('HNO4',   'swri',   'SWRI/HNO4.txt'),
    ('N2O5',   'swri',   'SWRI/N2O5.txt'),
    ('NO3',    'swri',   'SWRI/NO3.txt'),
    ('CH2N2',  'mpi',    'MPI/CH2N2_McMillan(1966)_298K_259.0-469.8nm.txt'),
    ('H2SO4',  'mpi',    "MPI/H2SO4_Farahani(2019)_0K_118.09-341.44nm(calc).txt"),
    ('S4',     'mpi',    'MPI/S4_BillmersSmith(1991)_723-843K_425-575nm.txt'),
    ('C2H3CN', 'mpi',    'MPI/C2H3CN_Eden(2003)_298K_113-320nm.txt'),
    ('C2H5SH', 'mpi',    'MPI/C2H5SH_ClarkSimpson(1965)_298K_162.0-249.0nm.txt'),
    ('CH3SOCH3', 'mpi',  'MPI/CH3SOCH3_Drage(2006)_298K_115-340nm.txt'),
    ('C3H4',   'mpi',    'MPI/C3H4_Nee(2008)_298K_105-220nm.txt'),
    ('C4H8',   'mpi',    'MPI/c-C4H8_RaymondaSimpson(1967)_298K_107.6-174.8nm.txt'),
    ('C4H6',   'mpi',    'MPI/C4H6-1_NakayamaWatanabe(1964)_298K_106.2-167.2nm.txt'),
]


def numeric_rows(path):
    """Yield each line's floats, for lines where every token parses as a float."""
    with open(path, errors='replace') as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                vals = [float(p) for p in parts]
            except ValueError:
                continue
            yield vals


def read_leiden(path):
    rows = [v for v in numeric_rows(path) if len(v) in (3, 4)]
    arr = np.asarray(rows)
    wn = 1e7 / arr[:, 0]                 # nm -> cm-1
    cs = arr[:, 1:].max(axis=1)          # max over absorption/dissociation/ionization
    data = np.column_stack([wn, cs])
    return data[np.argsort(data[:, 0])]


def read_swri_total(path):
    """First contiguous numeric block only (the 'Total'/sole first channel
    table) -- SWRI files have further branching-ratio tables later in the
    same file that must NOT be swept in."""
    rows = []
    started = False
    ncols = None
    with open(path, errors='replace') as fh:
        for line in fh:
            parts = line.split()
            vals = None
            if len(parts) >= 2:
                try:
                    vals = [float(p) for p in parts]
                except ValueError:
                    vals = None
            if vals is None:
                if started:
                    break
                continue
            if not started:
                started = True
                ncols = len(vals)
            if len(vals) != ncols:
                break
            rows.append(vals)
    arr = np.asarray(rows)
    wn = 1e8 / arr[:, 0]                 # Angstrom -> cm-1
    cs = arr[:, 1]                       # first value column = Total (verified per-file)
    data = np.column_stack([wn, cs])
    return data[np.argsort(data[:, 0])]


def read_mpi(path):
    rows = list(numeric_rows(path))
    arr = np.asarray(rows)
    wn = 1e7 / arr[:, 0]                 # nm -> cm-1
    cs = arr[:, 1]
    data = np.column_stack([wn, cs])
    return data[np.argsort(data[:, 0])]


READERS = {'leiden': read_leiden, 'swri': read_swri_total, 'mpi': read_mpi}

wn_grid = 1e8 / np.loadtxt(WN_GRID_FILE)   # AA -> cm-1 (nwreal_7949 is in Angstrom)
wn_grid = np.sort(wn_grid)                 # sort AFTER converting -- see module docstring


def fmt_tmol_lines(tmol):
    lines = []
    for i in range(0, len(tmol), 4):
        chunk = tmol[i:i + 4]
        lines.append(' '.join(f'{t:.1f} ,' for t in chunk))
    return lines


for molid, kind, relpath in SPECIES:
    path = f'{CROSSSECS}/{relpath}'
    data = READERS[kind](path)
    # de-duplicate wavenumbers, keep max cross section (same rule as production script)
    order = np.lexsort((-data[:, 1], data[:, 0]))
    data = data[order]
    _, uniq_idx = np.unique(data[:, 0], return_index=True)
    data = data[uniq_idx]

    molmass = Formula(molid).mass * u.u
    fact = (u.cm**2 / molmass).to(u.cm**2 / u.g).value
    cs_g = data[:, 1] * fact

    interp = PchipInterpolator(data[:, 0], cs_g, extrapolate=False)
    cs_interp = interp(wn_grid)
    cs_interp = np.where(np.isfinite(cs_interp), cs_interp, 0.0)
    cs_interp = np.maximum(cs_interp, 0.0)

    ktemp = len(TMOL)
    nwnos = len(wn_grid)
    out_path = f'{OUTDIR}/{molid}_UVonly_7949.dat'

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
        for wn, cs in zip(wn_grid, cs_interp):
            if wn == 0.0:
                continue
            cs_str = ' '.join(f'{cs:.3e}' for _ in range(ktemp))
            fh.write(f'{wn:.3f} {cs_str}\n')

    nonzero = (cs_interp > 0).sum()
    print(f'{molid}: {nwnos} pts written to {out_path}, '
          f'{nonzero} nonzero on MARCS grid, '
          f'source range {data[:,0].min():.1f}-{data[:,0].max():.1f} cm-1')
