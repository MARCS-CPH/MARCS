# Planet surface boundary: temperature clamp (SURFAC=1) vs energy-conserving surface (2026-11-23)

Handover for the agent that picks up the surface-boundary work. Nothing has been implemented yet: this
document is the exploration plan. Out of scope: the coupled chemistry solver (`toy_mixing/NONEQ_REFACTOR_STATUS.md`),
the water inventories (`WATER_INVENTORY_HANDOVER.md`), the archean chemistry stall (job 768509 diagnostic, see the status doc).
Line numbers refer to `marcs_noneq.f` (the candidate future `marcs.f`); the same code is expected in the tracked
`marcs.f` at an offset of roughly -6 to -17 lines, verify before patching.

## 1. The user's goal
Today `SURFAC=1` clamps the temperature of the last layer to `TSURF`. In the user's experience clamping one layer to a
fixed value hurts convergence. The idea: keep the last layer's temperature free and add the energy the surface
receives to that layer as an extra flux/heating, analogous to the species surface fluxes in the mixing solve, or by
adding a zeroth / NTAU+1 (surface) layer. Wording used: "surface_flux*albedo into the last layer (the layer just
above the surface)".

Open interpretation questions (ask the user before coding):
1. Which flux heats the layer? Physically the reflected part, albedo x F_down, is already in the radiative transfer
   (section 2); the part that is missing is the absorbed part, (1 - albedo) x F_down. Our proposal: heat with
   (1 - albedo) x F_down. Not yet confirmed by the user.
2. Does the surface also emit sigma*T^4 upward (T from the layer or a surface value), or is heating alone a first version?
   Our proposal: heating only first.
3. Change the radiative lower boundary (J = B(T) "deep diffusion", section 2) or leave it?

## 2. What the current mechanism does (read from the code, not yet tested)
Input: `SURFAC= TSURF= ALBEDO= NCH4=` is read at line ~6038 (format 1235) into `irrinp, tbottom, reflect, nch4`.
`irrinp>0` is the clamp. Naming trap: `HSURF` / `Y1` / the `TSURF=TF(TAU(1))` at line ~7871 refer to the TOP
boundary (stellar surface), not the planet surface; the planet surface temperature is the variable `tbottom`.

| What | Where (marcs_noneq.f) | Behaviour |
|---|---|---|
| Temperature equation, bottom layer | ~9239-9241 | `RT(NTAU) = -(TT(NTAU)-tbottom)`, `TTT(NTAU,NTAU)=1`: the radiative-equilibrium row of the last layer is REPLACED by the clamp. Only reached if `K <= MIHAL` (`IF (K.GT.MIHAL) GO TO 145` above it). |
| Re-assert before MATINV | ~9879-9915 | Gauss-Jordan elimination of column NTAU with the clamp row as pivot, then row NTAU zeroed. Runs for any MIHAL. The comment records an earlier half-applied version giving a -14651 K correction on iteration 1. |
| Albedo in the RT | ~10729-10732 (SCATTR) | `Pstar(JTAU) += reflect*Pstar(JTAU)`, then clipped to `[0, PSTARBND]` (the TOP-boundary beam strength). So the reflected beam is a local boost at one layer, and the clip can cut it. |
| Lower RT boundary | ~9206-9215 (`LOWER BOUNDARY`) | `XJ2=-1`, `J = B(T)` at the bottom, i.e. "an optically thick atmosphere continues below". This already acts as an unlimited heat source from below. |
| Absorbed stellar flux | nowhere | (1-albedo) x F_down at the bottom is never handed to any layer; the clamp absorbs whatever imbalance arises. |

Who uses it: 75 decks have `SURFAC=1` (all K2 `*_ocean_*` decks, the archean series, the carbon-burial GOE paper
runs, `full_network_archean_MIHALtest_*`). Of the benchmark/developing decks only `earth_archean` has `SURFAC=1`;
modern Earth, K2 planet, K2 virga and the developing decks have `SURFAC=0` (no surface at all; `TSURF` and `ALBEDO`
are read but unused). Therefore: do NOT change `SURFAC=1` behaviour. Add a new value (`SURFAC=2`).

Existing evidence that the clamp is a problem (found, not re-checked): the comment in
`benchmark_models/k2_18b_planet/slurm_script` says SURFAC=1 "was found to destabilize the adjacent layer under the
current SOLVE code, while SURFAC=0 settles cleanly" (test `K2_18b_Irr_spec_ocean_eq_surfac0_verify`, 2026-08-29).

## 3. Working hypotheses (all UNVERIFIED, each with a test in section 5)
- H1: the clamp is a hidden heat reservoir; the column can end with outgoing flux != absorbed stellar flux.
- H2: results depend on MIHAL because the anchor is only applied for `K<=MIHAL` while the elimination block always runs.
- H3: the albedo clip at `PSTARBND` removes part of the reflection at some wavelengths.
- H4: the clamp poorly conditions the neighbouring layer (documented instability, convergence trouble).
- Lead only: the printed `Current Teff...` is about 325 K (archean, SURFAC=1), about 330 K (modern Earth, SURFAC=0),
  about 308 K (K2 planet), versus equilibrium temperatures of 239-255 K for the Earth/Sun cases. The normalisation of
  that quantity (`TEFFP=TEFF*(FTOT/FLUX/PI)**.25`, line ~9355) was NOT checked, so this proves nothing yet.

## 4. Proposed design (needs user approval of the interpretation questions)
`SURFAC=2`: keep the radiative-equilibrium (and convective) row of layer NTAU, remove the clamp and the elimination
block, and add heating `q = F_abs * g / (P_NTAU - P_NTAU-1)` [erg g-1 s-1] with `F_abs = (1-albedo) * F_down(NTAU)`,
`F_down` being the downward stellar flux at the bottom (from `Pstar`; check units and the `f_irrad` dilution). The
temperature of the bottom layer then emerges. `TSURF` becomes unused (or only the starting guess). Later, if the
energy audit demands it: replace the `J=B` boundary by an upward intensity `(1-A)*B(T_s)` + reflected part, with `T_s`
from a surface energy balance (the "extra layer" version). Convection (mixing length) must carry the heat up from a
thin heated bottom layer; check the bottom layer becomes convective.

## 5. Suggested runs
Fast cases first. `NONEQ=0` decks need no chemistry solve (2-6 minutes): `input_files/K2_18b_Irr_spec_ocean_eq.input`
(SURFAC=1, TSURF=300, ALBEDO=0.07, MIHAL=72, ITEMAX=5; script `execute_scripts/runmarcs_K2_18b_Irr_spec_ocean_eq`,
plus `..._surfac0_verify`) and the archean deck with `NONEQ=0` (also tests the tuned elabund,
`data/elabund_earth_archean_tuned*.dat`). Generator for test runs: `toy_mixing/noneq_rename/make_bm_runs.py`
(`ONLY=<tag prefix>`; each run needs its own `out_dir`).

| ID | Run | Question answered |
|---|---|---|
| A0 | Add a diagnostic print (no physics change): incoming stellar flux, outgoing TOA flux (thermal+reflected), net flux through the bottom boundary, per iteration. Run SURFAC=0 and 1 on the K2 ocean eq deck and on archean NONEQ=0. | H1: is energy conserved with the clamp? |
| A1 | SURFAC=1, same restart, `TSURF` = 250 / 288 / 320 K | H1: if outgoing flux changes with TSURF while the stellar input is fixed, the clamp is a reservoir |
| A2 | SURFAC=1 with `MIHAL` = NTAU vs 10 vs 1 (decks `full_network_archean_MIHALtest_*` exist from earlier work; results unknown to us) | H2 |
| A3 | SURFAC=1 with `ALBEDO` = 0 / 0.3 / 0.9; print `Pstar(JTAU)` before/after the clip and count clipped wavelengths | H3 |
| A4 | Iterations to converge and the largest per-iteration T correction, SURFAC=0 vs 1, K2 ocean eq | H4 |
| B1 | After implementing: SURFAC=2 vs 1 on A0's cases: does it converge faster, what T_bottom emerges? | the main test |
| B2 | SURFAC=2 modern-Earth-like (1 AU Sun, albedo 0.30): is T_bottom near 288 K without being told? | energy balance sanity |
| B3 | SURFAC=2 scans: albedo 0.07 / 0.3 / 0.5; insolation (SPECSC 0.82 / 1.0); CO2 (1% vs 10%) | T_bottom responds in the right direction (the point of the change) |
| B4 | SURFAC=2 on the archean tuned deck (pre-GOE, SURFAC=1 today) and one K2 ocean case | regressions on the decks that use a surface |

## 6. Conventions (standing rules)
- Never push (SSH needs a passphrase; the user pushes). Commit only when asked, never `git add -A`.
- Edits to `marcs.f`: `marcsh lint` / `marcsh build` (`/groups/astro/tbalduin/marcs_harness/marcsh`), fixed form,
  column 72, declare every variable, comments at most 2 lines, show the proposed diff and wait for approval.
  `marcsh lint --base` needs a baseline file whose name ends in `.f`.
- Never leave `marcs.f` or `parameter.inc` modified after a build; builds with other settings (e.g. `NWREAL`) go from a
  scratch copy with `-I`.
- Every run gets a new `out_dir`; do not write into `benchmark_models/`. The Fortran log is buffered when it goes to a
  file: set `FORT_BUFFERED=false` in the job script when you need progress.
- node682 has a full `/dev/shm`; submit with `--exclude=node682`. `astro2_short` can queue for hours: `astro2_long` often starts sooner.
- Run scripts exit 0 even if MARCS stopped early: read the log (`marcsh check` where it applies).
- The Newton guard (`toy_mixing/guard/`) now stops a run after 5 consecutive rejected corrections (exit code 3).

## 7. Related open items (not part of this task)
Archean chemistry stall at chunk 10 (job 768509 queued), NaN T correction with the coupled chemistry on the CNDNSA+VIRGA
developing decks, stored benchmark references not reproducing, committing the refactor. See `toy_mixing/NONEQ_REFACTOR_STATUS.md`.
