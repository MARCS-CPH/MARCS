# Three water inventories in the CNDNSA / LATENT / moist-convection chain (2026-11-20)

Handover for the agent that picks up the moist-convection / water-consistency work. Read
`MOIST_CONVECTION_STATUS.md` first (FFL/FFC double count, moist-Cp code in `marcs_convec_test.f`).
Out of scope: Kzz (handled elsewhere), the clouds/VIRGA-opacity track, and the noneq solver
internals (status: `toy_mixing/NONEQ_REFACTOR_STATUS.md`).

## The problem
Inside one MARCS outer iteration, three different water inventories are used and nothing makes
them agree:

| # | Inventory | Where it lives | Who consumes it |
|---|---|---|---|
| 1 | Kinetic (noneq) water: whatever chemistry + Kzz mixing leave after the solve | written back to `ppallmol` by `noneq_return` (`marcs_noneq.f` ~15178; in `marcs_convec_test.f` it is `krome_return`, ~15476-15484) | opacities; the moist **saturation test** `h2o_pfrac = ppallmol(k,ih2o)/psum` (`marcs_convec_test.f:8636`) |
| 2 | GGchem *equilibrium* water at (T,p), with CNDNSA-depleted elements | `call ggchem(k, T, p, .False.)` inside `calc_adiaindex` (`marcs_convec_test.f:8499`) | the dry `cpR`, `gamma`, `xmmw` that the moist term is added to |
| 3 | VIRGA sedimentation flux `sedflux_in(k,i)` | read in OSTABLOOK (`marcs_noneq.f:13415`) | latent flux `FFL(k)` (SOLVE, `marcs_noneq.f` ~9440-9455; `marcs_convec_test.f:9500`) |

CNDNSA element depletion (`abundforgg_depleted`, `marcs_noneq.f:13551-13554`) only feeds GGchem
(inventory 2), never the kinetic solve. So the saturation test sees #1, Cp sees #2, the latent
heat sees #3.

## Evidence (own analysis, 2026-11-20; reproduce with `python_scripts/h2o_saturation_check.py <output_dir>`)
Coupled-scheme Earth runs (70 layers). S = p_H2O/p_sat(T), Buck 1996, ice below 273.15 K.

| Run (`toy_mixing/runs/`) | H2O vmr (min = max) | layers with S>1 | max S |
|---|---|---|---|
| `owndls_noflux_xsec3` (no CNDNSA) | 2.5e-3 | 8 | 1.6 |
| `owndls_cndnsa_latent` (CNDNSA+LATENT, 20 it) | 2.0e-3 | 14 | 2.6 |
| `owndls_cndnsa_latent2` | 2.2e-3 | 9 | 4.9 |
| `owndls_dryH2O` (H2O pinned 4e-6 aloft) | 4.0e-6..4.4e-6 | 0 | 0.003 |
| legacy-scheme `developing_models/..._cndnsa_latent_moistconvec` (13 Nov) | 4.6e-4..2.6e-3 | 3 | 1.4 |

- The profile is exactly uniform in the coupled runs: mixing time H^2/K is 0.3 yr (0.1 bar) to
  1.4 yr (1 bar) for the 540-yr solve, there is no condensation term in the kinetic ODE and the
  boundaries are closed. Whatever CNDNSA did to the start state is erased.
- The 13 Nov moist-convection test was made with the legacy per-layer scheme, which mixed weakly
  and kept some structure. The coupled scheme is now the only one: that test must be repeated.
- "Saturated" layers (588 moist engagements in layers 64-69) coincide with the supersaturated
  band of the uniform profile, so part of the moist effect is likely an artifact.
- A saturation ceiling alone cannot give a dry stratosphere: the minimum saturation mixing
  ratio along our columns is 0.8e-3..1.7e-3 (coldest troposphere 240-250 K at 0.3-0.45 bar),
  versus 3e-6..7e-6 for the real tropical cold point (190-195 K at 0.1 bar).
- Likely link to the ~10x ozone deficit: H2O 600x too high aloft is the main O3 sink (the
  H2O-pin test raised the O3 peak from 5.5e-7 to 5e-6).
- Unexplained, not traced: `cndnsa_diag.dat` (`..._cndnsa_latent/output`) shows `f_gas` for H =
  0.924 in layers 1-56 (unsaturated, 214-260 K) and 0.28-0.29 at layers 61-66. A uniform 7.6%
  depletion in unsaturated layers looks like a bookkeeping artifact of the VIRGA coupling.

## What VIRGA can and cannot do here
VIRGA (python `virga` 2.0.1) is an Ackerman-Marley-type sedimentation model. Its H2O vapour
pressure (`pvaps.py:H2O`) branches to ice below 273.16 K, so it can saturate and sediment ice
where the column is cold enough. It only sees the model T, P and Kzz and returns element
depletion for GGchem: it is a diagnosis on top of the model temperature structure, not a
Lagrangian cold trap, and the model has no cold tropopause (see above).

## Decisions made with the user
- End goal: self-consistent water. Evaporation and photolysis of water stay in the photochemical
  solver; no hand-tuned profiles in production.
- First test case (agreed as a good start): a **predetermined water profile with water removed
  from the kinetic evolution**, so inventories 1, 2 and 3 all see the same water.
  - Not a literally waterless network: removing the H2O reactions would delete HOx production
    (O(1D)+H2O, H2O photolysis). Keep H2O as a reactant at its prescribed abundance and freeze
    it, as is already done for the bath species M (`idx_M` in `noneq_dlsode_*`: ydot = 0, not
    mixed, excluded from the closure sum; generalise to a frozen-species mask).
  - INITAB as it is cannot express a vertical profile (one mixing ratio per species, last line
    wins, 3rd column = bottom N layers) and only resets at the start of each outer iteration.
    The profile needs a small new input (e.g. P-vs-vmr table, log-P interpolation) or an INITAB
    extension.
- Kzz: leave alone.

## Open questions
1. Which water should the Cp/GGchem call in `calc_adiaindex` use once water is prescribed?
2. How should `FFL` be computed with a prescribed profile (condensation rate implied by the
   profile, or drop FFL where the profile is prescribed)?
3. Repeat the moist-convection + gating test (`developing_models/..._moistconvec_gated`, not
   examined) under the coupled scheme.
4. Trace the uniform `f_gas(H)=0.924` depletion.

## Conventions
No pushes (the user pushes from their own terminal); never leave `marcs.f` or `parameter.inc`
modified after a build (back up, build, verify restore); code comments at most 2 lines; plots to
`/groups/astro/tbalduin/plots/`; Earth builds: use `krome/build_earth_jacfix` objects but the
current cross-section files (see the pitfalls in `toy_mixing/NONEQ_REFACTOR_STATUS.md`); do not run
three Earth jobs on one node (OOM).
