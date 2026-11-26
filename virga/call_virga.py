#!/usr/bin/env python3
"""Entry point marcs.f invokes via system('python3 call_virga.py').

Runs virga against the marcs2virga.dat / marcs_wnos.dat files marcs.f just
wrote in the current working directory (the SLURM job's scratch dir), and
writes virga2marcs.dat back into the current working directory for
marcs.f's OSTABLOOK to read.
"""
import os
import sys
import importlib.util

# needed for `import virga.justdoit` (the ext_mmr-patched copy lives
# here; some other, unpatched virga install must sit earlier on
# sys.path by default, since appending this instead of inserting it
# at position 0 made jdi.compute() reject ext_mmr -- confirmed
# 2026-11-26). Insert at 0 so this repo's virga.justdoit always wins.
VIRGA_MARCS_REPO = '/groups/astro/tbalduin/virga_marcs'
sys.path.insert(0, VIRGA_MARCS_REPO)

# Load virga2marcs_script.py by its explicit path next to this file,
# not via `import` -- VIRGA_MARCS_REPO has its own same-named copy,
# which a bare import silently preferred over this repo's (confirmed
# 2026-11-26, cost an entire fsed sweep: three runs all silently used
# the old copy). This is immune to sys.path order.
_here = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    'virga2marcs_script', os.path.join(_here, 'virga2marcs_script.py'))
_v2m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_v2m)
run_virga = _v2m.run_virga

# precomputed Mie scattering coefficients (784MB, kept outside git; see virga/Mie_data)
MIE_DATA_DIR = '/groups/astro/tbalduin/virga_mie_data'

if __name__ == '__main__':
    run_virga(input_dir='.',
              mieff_dir=f'{VIRGA_MARCS_REPO}/virga/',
              mie_data_dir=MIE_DATA_DIR)
