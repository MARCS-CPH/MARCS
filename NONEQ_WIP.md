# Branch `noneq-coupled-wip`: coupled chemistry + vertical mixing (work in progress, 2026-11-23)

This branch replaces the per-layer KROME solve + explicit mixing in `marcs.f` by ONE stiff ODE for the whole
column (DLSODE, banded analytic Jacobian from KROME, mixing inside the ODE). The old per-layer solver stays on `main`
(and tag `pre-noneq`). Status doc with details and run names: `toy_mixing/NONEQ_REFACTOR_STATUS.md`.

## What works (tested)
- Coupled solver on Earth (94 species), K2-18b (67 species) and a network without the third-body species M.
  K2 regression: 66/67 species agree with the earlier run to <=1.6e-4 (only a species with no sink differs).
  No-M test: max relative difference 1.5e-4.
- `NONEQ=0` benchmarks (K2 planet, K2 virga, K2 star, solar, Trappist): identical to the old code.
- Earth benchmark decks (`itmax5`, `mihal10`) and the 7 `developing_models` decks run to completion
  (the developing decks need a KROME build with 16937 photo bins and `NWREAL=16937`).
- INITAB: per-layer profiles (3rd column = first layer counted from the top), optional 4th column = floor fraction.
- Surface-flux file `noneq_surface_flux.dat` with an optional layer column; source/sink audit in the log.
- Newton-correction guard v2: NaN/Inf and relative-size test, `GUARD` log line, stop after 5 rejected corrections in a row
  (exit code 3). Fixes the Trappist benchmark that never updated its model.
- Settings in `noneq_parameter.inc` (compile time); the deck's `tMAX`/`DTMIN`/`DTMAX` are now end time / first step / step cap.

## Known problems / not done
- Archean (133 species) at a realistic `tMAX`: unstable from t ~ 5e5 s (peroxy radicals in the top layers go negative);
  under investigation.
- Four CNDNSA=1 + VIRGA=1 developing decks: NaN temperature correction at layer 1 on iteration 2 with the coupled chemistry
  (old code: none). The guard rejects it and the run continues.
- Removed relative to `main`: AQUA, BURIED, TOYMIX, two-cycle mixing, the per-layer path. Decks with `CNVCRI=0`
  now run the coupled solve; stored benchmark outputs of those decks are not regenerated.
- Stored reference outputs of several benchmarks do not reproduce even with the old code (K2 star 14 K, solar 51 K, Trappist 50 K).
- Not committed here: moist convection / water inventories (see `WATER_INVENTORY_HANDOVER.md`), surface-boundary rewrite
  (see `SURFACE_BOUNDARY_HANDOVER.md`), freezing of species.

## Building
1. Regenerate the KROME network with the patched generator (`krome/kromeobj.py`; the generated `krome/build_*` dirs are not tracked).
   The earlier builds contain a wrong Jacobian (duplicate `elseif`, uninitialised `pdj`).
2. `make KROME_BUILD=krome/build_xxx EXEC=marcs_xxx` (the Makefile now rebuilds when `parameter.inc` or `noneq_parameter.inc` change).
3. New input/output files are named `noneq_*` (old names are still read as a fallback).
