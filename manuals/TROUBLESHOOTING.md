# Running a MARCS model, and what goes wrong

For first-time setup (compilers, GGchem/krome build steps, required
data files) see `../README.md`. For the full parameter reference see
`parameters_list.txt`. For physics-module explanations see
`PHYSICS_OVERVIEW.md`. This document is about actually running a model
day-to-day and diagnosing failures — written from a session
(2026-08-20/21) that hit essentially every failure mode below firsthand
while doing a large `marcs.f` cleanup pass.

## How a run actually works

Every `execute_scripts/runmarcs_*` script follows the same shape:

1. Two directory variables near the top: `in_dir` (where the restart
   archive/reference data this run continues from lives) and `out_dir`
   (where this run's results get written). **These being the same
   directory, or `out_dir` pointing at anything git-tracked, is the
   single most common way to lose data** — see "Overwriting tracked
   output" below.
2. `input_parameters=$CWD${input_folder}<deck>.input` — the actual
   `.input` deck used, pulled from `input_files/` (or a subdirectory of
   it).
3. Everything gets staged into `$TMPDIR` on `/dev/shm` (fast local
   scratch on the compute node), the compiled `marcs` binary runs there,
   and results get copied back to `out_dir` afterward.
4. Submit with `sbatch execute_scripts/runmarcs_<name>` **from the repo
   root** — every path inside is relative to wherever the script is
   submitted from (`CWD=$(pwd)`), not to where the script file itself
   lives.

### The 4 canonical benchmark models

`solar_type_star`, `k2_18b_star`, `k2_18b_planet`, `trappist1_star` each
live in `benchmark_models/<name>/` — this is git-tracked *reference*
data (restart archive, original output). Each has its own
`benchmark_models/<name>/output/` subdirectory (gitignored) that a run
writes fresh results into. **Never point a run's `out_dir` back at
`benchmark_models/<name>/` itself** — only `benchmark_models/<name>/output/`.
This split exists specifically because every one of these 4 launch
scripts used to write output on top of (or into the same directory as)
the tracked reference data, and it happened for real, more than once, in
the session that wrote this document.

### Running your own model

Copy an existing script that resembles what you want (similar
`NONEQ`/`PHOTO`/`IRRAD` settings), point `input_parameters` at your
deck, and pick a fresh `out_dir` that isn't a path anyone else's script
also uses — e.g. `models/my_test_name/`, never `benchmark_models/*/`
directly.

## Common errors and how to fix them

### 1. Deck is out of date for the current `marcs.f` (the most common problem)

`marcs.f` has had a lot of dead input parameters removed over time
(`IDUST`, `sampling`/`OS-ODF`, `FACPLY`, `METBL`, `MOLTSUJI`, `ISPEC`,
`lops`, `KROMEO`, `JUMP`, `NLTE`, `ILINE`, `sunz`, `abundatms_inp` — see
`parameters_list.txt` for the full removed-parameter list with exact
column positions). A `.input` deck that predates one of these removals
no longer lines up with what the code expects. **This shows up in two
very different ways depending on where the stale field is:**

**(a) Fixed-column fields — silent misalignment, not a crash.** Fortran
reads these by character position, not by parsing "name=value" text. If
a deck still has an extra field the code no longer reads, every field
after it on that line gets read from the wrong columns — you get a
plausible-looking but wrong number, or a downstream crash somewhere
unrelated once that garbage value gets used (e.g. a nonsensical `METPE`
value). **No error message points you at the real cause.**

**(b) NAMELIST fields (`&abundances`, `&inputos`, `&outlist`) — a hard
crash with a specific, recognizable error:**
```
forrtl: severe (19): invalid reference to variable in NAMELIST input,
unit 5, file /dev/shm/tmp.XXXXXX/mxms7.input, line 24, position 19
```
This means the deck's `&abundances` (or other namelist) line has a
field name (e.g. `sunz=`, `abundatms_inp=`) that the current `marcs.f`
no longer declares in that namelist. Unlike fixed-column fields,
NAMELIST reads are normally order-independent and would silently accept
missing fields — but ifort's reader is strict about *unrecognized*
field names and aborts immediately. This is actually the easier case to
diagnose, precisely because it fails loudly.

**How to fix**: compare your deck's header lines against a known-current
one — `input_files/marcs_sun.input` in this repo is kept up to date and
is a good reference. For example, the current (correct) first two lines
and `&abundances` line look like:
```
ITEMAX=150     NPRINT=  0     NEWMOD=  0     NOARCH=  0     JONTYP=3
MOLOLD=  1     PP/SPH=PP
 &abundances  zscale=    1.00000
```
An out-of-date deck instead has (note the extra `IDUST`, `JUMP`, and
`OS/ODF` fields, and `sunz=`/`abundatms_inp=` in the namelist):
```
ITEMAX=150     NPRINT=  0     NEWMOD=  0     NOARCH=  0     JONTYP=3       IDUST =0
JUMP  =  4     OS/ODF=OS      MOLOLD=  1     NLTE  =  0     PP/SPH=PP
 &abundances  sunz='y', zscale=    1.00000, abundatms_inp=  0.    12.000000000000
```
Fix by deleting exactly the obsolete field(s) and nothing else — for
fixed-column fields this usually means removing a fixed number of
characters (check `parameters_list.txt` for the exact field width) and
letting the remaining columns shift left; for namelist fields, just
delete `name=value,` and its neighboring comma. Always verify the
result against `parameters_list.txt`'s current format before running,
not just by eye — column-count mistakes are easy to make and produce
exactly the "no error, wrong answer" failure this section opened with.

### 2. Overwriting tracked benchmark output

Symptom: `git status` (or `git diff`) shows a tracked file like
`benchmark_models/solar_type_star/marcs_output.arciv` as modified after
you just wanted to test something. This means a script's `out_dir`
pointed at tracked reference data. Undo with `git checkout --
<path>` (only if you haven't committed — check first), then fix the
script to write somewhere untracked before rerunning. See "The 4
canonical benchmark models" above for the intended pattern.

### 3. GGchem template not found

`marcs.f`'s `GGCHEM` routine opens `marcs2ggchem_template.in` with no
path prefix, so it must be staged at `$TMPDIR`'s top level specifically
— being present only inside the copied `GGchem/` subdirectory isn't
enough. If GGchem calls are failing with a file-not-found-style error,
check the run script has a line like
`cp GGchem/marcs2ggchem_template.in $TMPDIR/` (not just
`cp -r GGchem $TMPDIR/`).

### 4. VIRGA staging incomplete

If running with `VIRGA=1`, note this coupling has had staging gaps
before (a stale input deck and a missing template file, unrelated to
VIRGA's own physics) — if a VIRGA-enabled run fails early, check first
whether `call_virga.py` and its dependencies actually got copied into
`$TMPDIR`, before assuming the physics is broken.

### 5. Stale cached `krome_flux_rad.dat`

If a run script copies in a `krome_flux_rad.dat` left over from a
*different* model's grid, you'll see an error like:
```
ERROR: krome_flux_rad.dat ended before supplying all ntau*nwreal
records needed for this run's grid...
This file was carried over from a previous run with a different
LOGTAU/wavelength grid and no longer matches
```
Fix: comment out the line copying `krome_flux_rad.dat` into `$TMPDIR`
in the run script so it regenerates fresh instead of reusing the stale,
mismatched cache.

### 6. Missing or incompatible restart archive

Every run — even a "fresh" one (`NEWMOD=1`) — reads a restart archive
(`arcivaaa.dat`, staged from `${in_dir}<name>.arciv`) at startup via
`SUBROUTINE OLDSTA`. If this file is missing, or was written by a
version of `marcs.f` with a different archive layout (e.g. an old
archive missing the `FOLD` block a `MOLOLD=1` deck expects to read), the
run will fail immediately with a low-level Fortran I/O error rather than
a descriptive message. Make sure `in_dir` actually points at a valid,
current-format `.arciv` file.

### 7. Run-to-run differences that are NOT a bug

This solver has inherent run-to-run floating-point nondeterminism
(confirmed empirically: running the exact same unmodified binary against
the exact same deck twice can produce Teff differing by roughly 0.2K,
likely from timing/thread-scheduling-dependent summation order). Don't
chase a ~0.1-0.2K discrepancy between two runs as a regression — it's
expected. Restarting every comparison from one fixed, never-touched
archive (rather than chaining from a previous run's slightly-noisy
output) makes comparisons effectively bit-exact in practice, and is the
right way to actually test whether a code change altered anything.

### 8. Genuine non-convergence (as opposed to just needing more iterations)

Watch the "Max corr. to T wanted" line in the run's log across
iterations. Ordinary slow convergence shows corrections shrinking
(unevenly, but trending down) toward `TCONV`. A **real, unconverged
oscillation** instead shows the same few correction values (often at the
same one or two layers) repeating in a fixed cycle without shrinking,
right up to `WARNING: Maximum number of iterations reached`. This is a
known, real failure mode in this code — a convective bistability bug
was found and fixed in `SOLVE` but, as of this writing, not yet ported
to `SOLVE_sph` (spherical-mode runs) — and at least one plane-parallel
case (`benchmark_models/k2_18b_star`) is currently exhibiting a similar
unresolved period-3 oscillation under investigation. If you see this
pattern, it's worth reporting rather than just re-running with a higher
`ITEMAX` and hoping.

### 9. Stale compiled binary

If `marcs.f` (or anything it depends on) was edited but `make` wasn't
rerun since, the run uses the old logic silently — there's no warning
that the binary predates the source. Check `ls -la marcs marcs.f` and
compare timestamps if a run behaves unexpectedly right after a code
change; rebuild with plain `make` if the binary is older.

### 10. `krome/MARCS_build/` not populated

`make` needs `krome/MARCS_build/` populated with `.o`/`.mod` files
before it can link `marcs`. If starting from a fresh checkout without
running KROME's own build step (`krome/compile_KROME.sh`, see
`../README.md`), `make` will fail at the link stage. This directory can
also be bootstrapped by copying it from another checkout that already
built the same network/version, if rebuilding from scratch isn't
practical right away.
