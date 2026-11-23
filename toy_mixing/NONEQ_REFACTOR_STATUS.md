# noneq refactor: status (2026-11-20)

Goal: make the coupled chemistry+mixing solve (one DLSODE ODE per column) the only scheme, rename
MARCS-owned `krome_*` names to `noneq_*`, and cut the run-time inputs. Not pushed; nothing in
production was changed (all work is in new files, every step verified against a reference run).

## Files
| File | What |
|---|---|
| `marcs_noneq.f` | working source (untracked). From `toy_mixing/noM_test/marcs_noM.f` by the scripts below |
| `noneq_parameter.inc` | compile-time numerics (rtol, atol, convergence test, chunk count, floor tau/band) |
| `toy_mixing/noneq_rename/step1a_rename.py` | identifiers, routines, module vars, log tags |
| `toy_mixing/noneq_rename/step1b_remove_legacy.py` | removes the per-layer path (see "Removed") |
| `toy_mixing/noneq_rename/step2_namelist_floor.py` + `noneq_dlsode_block.f` | namelist removal, floor, layer column; the block is the new solver source of truth |
| `toy_mixing/noneq_rename/check_fixed_form.py`, `wrap_fixed_form.py` | fixed-form column-72 check / re-wrap (ifort silently drops text past col 72) |
| `marcs_noneq.f.pre_step1b`, `.pre_step1d`, `.pre_step2` | backups before each step |
| `toy_mixing/noneq_rename/ref_step1b.f` | step-1b source, builds the namelist reference binary |
| binaries (MARCS dir) | `marcs_earth_noneq` (new), `marcs_earth_noneq_1b` (step-1b reference), built with `KB=krome/build_earth_jacfix toy_mixing/build_one.sh` |

## Behaviour changes
- `NONEQ=1` always runs the coupled solve. `MIXING=0` -> mixing scale 0. `CNVCRI<=0` -> no early stop.
- Deck: `tMAX` = end time, `DTMIN` = first DLSODE step, `DTMAX` = step cap (0 = none), `DTINC` ignored.
  Chunk ends are tMAX/2**(max_chunks-1), doubling. Warning if tMAX < 1e6 s (old decks).
- Removed: per-layer KROME loop, Thomas mixing, AQUA/BURIED transfer (networks with `*_BURIED`/`*_AQUA`
  species now stop with a message), TOYMIX, inner-loop passes, `inner_loop.nml`, hold/pin, prerelax,
  `noneq_full_output.dat`, `noneq_bins_rates.dat`. Those models stay on the older GitHub release.
- New files: outputs `noneq_final_output.dat` (Time column now real), `noneq_flux_rad.dat`;
  inputs `noneq_init_abund.dat`, `noneq_surface_flux.dat`, `noneq_flux_rad.dat` are read first, the old
  names are still accepted. Old slurm scripts that copy `krome_*_output.dat` will silently miss outputs.
- INITAB 4th column (0..1) = floor fraction of the listed mixing ratio during the solve. Smooth (C1)
  push starting at floor*(1+band), band 5% (`noneq_floor_band`): a held species sits ~1.05x the floor.
- Surface-flux file optional 4th column = layer per species (default bottom).
- INITAB profiles (2026-11-20): 3rd column = FIRST LAYER of the entry counted from the TOP (empty or <=1 = whole
  column; was 'layers from the bottom'). A species may have several lines; each layer takes the entry with the
  largest start layer at or above it (line order does not matter). Layers above the first entry stay unpinned.
  Startup log prints the applied profile (`INITAB H2O : L 1- 58= ...`). Tests: pushed draft prints the stepped
  profile, reversed line order gives identical output (max diff 0), single-line CH4 regression byte-identical.
  `krome_init_abund_earth_archean_co2pin5.dat` still uses the OLD bottom-N meaning (carbon-burial models, not
  supported here); the pushed layered-H2O draft was rewritten to starts 1/59/66 (backup `.bak_bottomN_convention`).
  INITAB only sets the START of each outer iteration: the long run flattens the stepped H2O to 2.26e-3 everywhere.
- INITAB otherwise only resets at the start of each outer iteration (CH4 pinned at 1.9e-6 ends at 1e-18..1e-28).
- `idx_M` = index of the third-body species M, 0 if the network has none (guards verified, see below).

## Verification (Earth, KB=build_earth_jacfix, reference-era cross sections)
| Step | Result |
|---|---|
| rename + M guards, outputs renamed | byte-identical to the reference run (runs `noneq_step1a`, `noneq_step1d`) |
| legacy path removed | byte-identical (`noneq_step1b`) |
| deck-driven settings vs namelist build, matched values, no flux and all fluxes | byte-identical states, 0/210 layers differ (`noneq_e0/e1_ref/new`), incl. after the floor smoothing |
| floor 0.8 on CH4 | 1.59e-6 everywhere vs floor 1.52e-6 (`noneq_f_floor`); sharp max(0,.) version collapsed DLSODE |
| layer column | CH4/SO2 source at layer 40 (`noneq_l_layer`) |
| no-M network (45 M reactions removed) vs same with a null M reaction | 1.5e-4 at rtol 1e-4, 8.7e-8 at rtol 1e-8 (solver roundoff) |
| prerelax A/B | ~25% fewer steps, end state within 2.6e-5: removed |
Pitfalls: `build_earth` lacks the Jacobian fix (use `build_earth_jacfix` objects) but the jacfix `.dat`
xsecs are stale (stage `toy_mixing/noM_test/xsec_as_of_ref/` for the 18 Nov reference, or `build_earth`
for current); three concurrent Earth jobs on one node get OOM-killed.

## Open
1. DONE 2026-11-20: K2 regression (run `noneq_final_k2`, ~25 min solver) vs B2: 66/67 species agree to <=1.6e-4; only CH3S (no sink, grows linearly) differs by 14% = ratio of end times (1e13 vs 8.6e12 s).
2. DONE 2026-11-20: no-M re-test with `marcs_noneq.f` (runs `noneq_final_noM`/`_nullM`): max rel diff 1.5e-4, 0/93 species >1%. Clean `make` with the committed Makefile also OK (Makefile still lacks the `noneq_parameter.inc` dependency).
3. Split the solver into its own file + small `marcs.f` patch; Makefile dependency on `noneq_parameter.inc`.
4. Update `NEW_SCHEME_OVERVIEW.md` (still describes the namelist and hold/pin) and the INITAB docs.
5. Commit the KROME Jacobian patches (`krome/kromeobj.py`, `krome/src/krome_user.f90`): hard dependency.
6. Benchmarks/developing decks with CNVCRI=0 (18) change meaning (now the coupled solve): regenerate.
7. User pushes from their own terminal (SSH passphrase).

## Git state (2026-11-20)
- Local commit `b45557b` (not pushed): INITAB per-layer profiles in the tracked `marcs.f` (+ rewritten layered-H2O draft).
  Tested with the legacy solver: 2-column files byte-identical old vs new (runs `legacy_l1/l2`, `legacy_l6/l7`),
  stepped H2O profile as expected and order independent (`legacy_l3/l5`; `legacy_l4` shows the old last-line-wins).
  Patch generator: `toy_mixing/legacy_initab/make_initab_patch.py`.
- Still uncommitted: `marcs_noneq.f`, `noneq_parameter.inc`, KROME Jacobian patches, README.md.

## Pre-push benchmark sweep (2026-11-23; runs toy_mixing/runs/bm_*, scripts toy_mixing/noneq_rename/make_bm_runs.py)
- NONEQ=0 benchmarks (planet, virga, k2star, solar, trappist): new code == old code (T identical; virga other columns <=1e-4). Stored reference outputs do not reproduce even with old code (planet 0.06 K, virga 0.6 K, k2star 14 K, solar 51 K, trappist 50 K).
- NONEQ=1: modern_earth itmax5/mihal10 and 7 developing decks run (dev needs krome/build_earth_wnos1_jacfix, 16937 photo bins, and NWREAL=16937 compiled from a scratch parameter.inc). 4 CNDNSA=1 dev decks: NaN T correction at layer 1, iteration 2, only with the coupled chemistry (old code: none).
- earth_archean (133 species, 13300 unknowns) with tMAX=3.2e10: stalled in chemistry chunk 10 of iteration 1, TIMEOUT after 4 h (min=-4e5, change 695 at chunk 9). Deck as-is (tMAX=1e3 = 1000 s solve) runs but is meaningless. UNRESOLVED.
- Decks with old tMAX semantics (archean tMAX=1e3, DTMAX=1e2) need new values; dev scripts exit 0 even when MARCS stops early (check the log).
- Newton guard v2: see toy_mixing/guard/ and memory project-marcs-newton-guard-v2.
