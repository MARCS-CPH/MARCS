# Physics overview (WIP)

This document is a work in progress and will be updated slowly over the
coming weeks. It collects explanations of MARCS's physics modules -- what
they do, which input parameters control them, and any current limitations
-- as a companion to ../README.md (which covers running the code) and
parameters_list.txt (the full input parameter reference, in this same
manuals/ directory).

### CIA (collision-induced absorption)
MARCS can include collision-induced absorption via the `LINCIA` input
parameter (0 = off, the default in every existing .input file; 1 = on).

When enabled, MARCS reads its list of active pairs from a fixed path,
`./data/cia_list.dat` -- this file is **not** created automatically, so you
need to put one there yourself before setting LINCIA=1. Format: first line
is the number of pairs, then one line per pair as
`SPECIES1  SPECIES2  /path/to/converted/pair/data/`.
`data/cia_pairs_all.dat` lists every CIA pair currently converted and
available (19 pairs) in exactly this format and can be used as a source to
copy the pairs you want into `data/cia_list.dat` -- it is a reference list
only and is not read directly by MARCS. `data/cia_pairs_all_with_ar.dat`
is the same but also includes the two Ar-bearing pairs (CO2-Ar, N2-Ar) --
only usable with an elabund file that tracks Ar as a species, otherwise
MARCS stops with "unrecognized species".

Pairs are matched generically by species name against GGchem's species
tables, so any pair present in the converted data can be used, not just the
originally-targeted N2/O2/CO2 set for Earth-type atmospheres.

### N2/O2 Rayleigh scattering
N2 and O2 Rayleigh scattering is included automatically whenever those
species are present in the atmosphere. There is no input switch for this --
it runs independently of LINCIA/CIA.

### conv_crit (chemistry integration length)
Part of the noneq input block (see parameters_list.txt). A value of exactly
`0` switches per-layer KROME chemistry integration to always run for the
full `tMAX`, instead of the default `min(vert_mix_time, tMAX)`. Any nonzero
value keeps the default behavior. This gives direct input-file control over
integration length independent of the vertical mixing timescale.

## Known physics gaps

### No water-vapor cold trap (missing condensation) -- VIRGA should fix this
The current non-eq columns hold H2O at an essentially **flat mixing ratio
across the whole depth grid**, regardless of local temperature -- e.g. a
modern-Earth test column showed ~2.54e-3 (2540 ppm) from a 230 K upper
layer down to the 288 K surface, with no falloff at all in between. Real
atmospheres don't do this: water vapor condenses out once it hits its local
saturation vapor pressure, which collapses sharply at low T (the "cold
trap") -- it's why Earth's real stratosphere is desiccated to a few ppm
despite a humid surface. There's currently no condensation physics coupled
into the non-eq column to enforce this, so GGchem/KROME just carries
whatever H2O partial pressure comes out of equilibrium chemistry at every
layer, uncapped by saturation.

This isn't a cosmetic issue -- it directly distorts the chemistry. **O(1D) +
H2O -> 2OH** is the dominant real-world tropospheric OH source, and it's a
fast, nearly gas-kinetic reaction. With H2O artificially high at every
altitude (including cold layers where it should be nearly absent), OH comes
out 10-75x too high in number density and 3-5 orders of magnitude too high
in mixing ratio through most of the column (found while investigating why
DMS never accumulates to a realistic level in the Earth-modern DMSP/DMS
tests -- see krome/networks/non_eq/full_network_DMS_v2_airsea.dat and
python_scripts/DMS_airsea_flux_calc.py). At the layers where OH is worst,
the OH + DMS abstraction channel alone destroys DMS 2-15x faster than it
can be produced -- a downstream symptom of this upstream H2O problem, and a
concrete illustration of how far this can propagate: anything that reacts
with OH (which is most of the network) inherits some of this distortion.

The in-progress VIRGA coupling (`virga/`, `call_virga.py`,
`virga2marcs_script.py`, invoked from marcs.f when `VIRGA=1`) is aimed at
exactly this gap -- `virga/virga_condensates.dat` already lists `H2O` as a
condensate species. Once wired in for the non-eq column (it's currently
used for cloud opacity, not yet feeding a depleted H2O profile back into
the KROME chemistry), it should let water actually condense out with
altitude/falling temperature instead of sitting at a fixed mixing ratio,
which should bring OH (and everything downstream of it) much closer to
realistic values.
