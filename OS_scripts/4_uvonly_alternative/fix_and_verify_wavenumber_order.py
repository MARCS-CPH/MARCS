#!/usr/bin/env python3
"""
Diagnostic + fix utility for the wavenumber-order bug (see the module
docstring in build_uvonly_flat.py for the full history). Two independent
uses:

1. verify_intermediate(path) -- check a *_UVonly_7949.dat intermediate text
   file's data rows are ascending in wavenumber.
2. reverse_intermediate(path) -- fix one by reversing its data rows in
   place (header block preserved). Safe because the wavenumber/cross-section
   pairing in these files is always correct -- interp(wn_grid) evaluates
   element-wise regardless of wn_grid's order, so a pure row-order bug is
   fixed by a pure row-order reversal, no recomputation needed.
3. verify_final_binary(path) -- check a converted output/combined/<X>/wn.dat
   Fortran-unformatted binary is ascending (this is what MARCS actually
   reads -- always verify the FINAL output, not just the intermediate, in
   case convert_opac itself ever reorders something).

Usage:
    # check one intermediate file
    python3 fix_and_verify_wavenumber_order.py verify-intermediate <path>

    # fix one intermediate file (reverses rows, keeps header)
    python3 fix_and_verify_wavenumber_order.py fix-intermediate <path>

    # check a final binary wn.dat (needs scipy)
    python3 fix_and_verify_wavenumber_order.py verify-binary <path>

    # batch-fix + report: edit SPECIES below to the affected list, run with
    # no arguments
    python3 fix_and_verify_wavenumber_order.py
"""
import sys
import numpy as np

INTERMEDIATE_DIR = '/groups/astro/tbalduin/marcs_opac_converter/combined_UV_IR/input_crossec_files/combined'
OUTPUT_DIR = '/groups/astro/tbalduin/marcs_opac_converter/combined_UV_IR/output/combined'

# Edit to the species you need to check/fix when reusing this script.
SPECIES = [
    'C2H', 'C2H3', 'C2H5', 'C3H3', 'CH2', 'CH3CHO', 'CH3NH2', 'CH3SH', 'HCO',
    'HNC', 'HNCO', 'NH2', 'HNO4', 'N2O5', 'NO3', 'CH2N2', 'H2SO4', 'S4',
    'C3H6', 'C4H4', 'C6H6', 'N2H4',
    'C2H3CN', 'C2H5SH', 'CH3SOCH3', 'C3H4', 'C4H8', 'C4H6',
]


def _split_header(lines):
    for i, l in enumerate(lines):
        if l.strip() == '/':
            return lines[:i + 1], lines[i + 1:]
    raise ValueError("no '/' header terminator found")


def verify_intermediate(path):
    with open(path) as f:
        lines = f.readlines()
    _, data = _split_header(lines)
    wns = np.array([float(l.split()[0]) for l in data])
    ascending = bool((np.diff(wns) > 0).all())
    print(f'{path}: {len(wns)} rows, {wns.min():.3f}->{wns.max():.3f}, '
          f'ascending={ascending}')
    return ascending


def reverse_intermediate(path):
    with open(path) as f:
        lines = f.readlines()
    header, data = _split_header(lines)
    first_wn = float(data[0].split()[0])
    last_wn = float(data[-1].split()[0])
    was_descending = first_wn > last_wn
    data_rev = list(reversed(data))
    with open(path, 'w') as f:
        f.writelines(header)
        f.writelines(data_rev)
    print(f'{path}: was {first_wn:.1f}->{last_wn:.1f} (descending={was_descending}), '
          f'now ascending, {len(data)} rows')


def verify_final_binary(path):
    from scipy.io import FortranFile
    ff = FortranFile(path, 'r')
    wn = ff.read_reals(dtype=np.float64)
    ff.close()
    ascending = bool((np.diff(wn) > 0).all())
    print(f'{path}: {len(wn)} pts, {wn.min():.3f}->{wn.max():.3f}, '
          f'ascending={ascending}')
    return ascending


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == 'verify-intermediate':
        ok = verify_intermediate(sys.argv[2])
        sys.exit(0 if ok else 1)
    elif len(sys.argv) == 3 and sys.argv[1] == 'fix-intermediate':
        reverse_intermediate(sys.argv[2])
    elif len(sys.argv) == 3 and sys.argv[1] == 'verify-binary':
        ok = verify_final_binary(sys.argv[2])
        sys.exit(0 if ok else 1)
    elif len(sys.argv) == 1:
        print(f'Checking + fixing {len(SPECIES)} intermediate files...')
        for sp in SPECIES:
            path = f'{INTERMEDIATE_DIR}/{sp}_UVonly_7949.dat'
            if not verify_intermediate(path):
                reverse_intermediate(path)
        print('\nDone. Now reconvert via 3_convert_opac/, then verify the final binaries:')
        all_ok = True
        for sp in SPECIES:
            path = f'{OUTPUT_DIR}/{sp}_UVonly/wn.dat'
            try:
                if not verify_final_binary(path):
                    all_ok = False
            except FileNotFoundError:
                print(f'{path}: not converted yet')
                all_ok = False
        print(f'\nAll ascending: {all_ok}')
    else:
        print(__doc__)
        sys.exit(1)
