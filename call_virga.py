#!/usr/bin/env python3
"""Entry point marcs.f invokes via system('python3 call_virga.py').

Runs virga against the marcs2virga.dat / marcs_wnos.dat files marcs.f just
wrote in the current working directory (the SLURM job's scratch dir), and
writes virga2marcs.dat back into the current working directory for
marcs.f's OSTABLOOK to read.
"""
import sys

VIRGA_MARCS_REPO = '/groups/astro/tbalduin/virga_marcs'
sys.path.insert(0, VIRGA_MARCS_REPO)

from virga2marcs_script import run_virga

if __name__ == '__main__':
    run_virga(input_dir='.',
              mieff_dir=f'{VIRGA_MARCS_REPO}/virga/',
              mie_data_dir=f'{VIRGA_MARCS_REPO}/Mie_data')
