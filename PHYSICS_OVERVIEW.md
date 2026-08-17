# Physics overview (WIP)

This document is a work in progress and will be updated slowly over the
coming weeks. It collects explanations of MARCS's physics modules -- what
they do, which input parameters control them, and any current limitations
-- as a companion to README.md (which covers running the code) and
parameters_list.txt (the full input parameter reference).

### CIA (collision-induced absorption)
MARCS can include collision-induced absorption via the `LINCIA` input
parameter (0 = off, the default in every existing .input file; 1 = on).

When enabled, MARCS reads its list of active pairs from a fixed path,
`./data/cia_list.dat` -- this file is **not** created automatically, so you
need to put one there yourself before setting LINCIA=1. Format: first line
is the number of pairs, then one line per pair as
`SPECIES1  SPECIES2  /path/to/converted/pair/data/`.
`data/cia_pairs_all.dat` lists every CIA pair currently converted and
available (21 pairs) in exactly this format and can be used as a source to
copy the pairs you want into `data/cia_list.dat` -- it is a reference list
only and is not read directly by MARCS.

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
