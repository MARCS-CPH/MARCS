# GGchem version comparison: vendored copy vs. current upstream

## Where each copy lives

| | Path | Bulk date | Compiled binary |
|---|---|---|---|
| Vendored (used by marcs) | `MARCS_CPH/MARCS/GGchem/src16` | 2025-03-12 | `ggchem16`, 2024-11-22 |
| Current upstream | `/groups/astro/tbalduin/GGchem/src16` | active git clone, `origin` = `github.com/pw31/GGchem.git`, branch `master` | latest local commit `f6c9f1b`, 2026-03-20 |

The upstream repo has no tags/releases/CHANGELOG; the only version markers are
git commit dates. Recent commit history relevant to the gap:
```
f6c9f1b use_SiO=F only takes effect in supersat, take some action in equil_cond after 15 unsuccessful iterations
05eda8f new strategy to avoid switching the same condensate on/off many times
2ca8950 remove H2SO4 tweak options from model_VenusStruc_Paul.in
af9b9ec Merge branch 'master' of https://github.com/pw31/GGchem
f0d2919 keep colours of condensates in Plot_structure.py
```

marcs.f calls GGchem as a separate OS process (`marcs.f:15589`,
`call system('./GGchem/ggchem marcs2ggchem.in > ggchem_out.txt')`), so the
"version in use" is entirely determined by whatever is compiled under
`MARCS_CPH/MARCS/GGchem/src16` — it does not automatically track upstream.

## File-by-file diff (`diff -rq` between the two `src16` trees, contents inspected directly)

**Unchanged — the call contract itself:**
- `ggchem.f` (`SUBROUTINE GGCHEM(nHges,Tg,eps,merk,verbose)`, the top-level
  gas-phase chemistry entry point) is **byte-identical** in both versions.
  If marcs.f were ever changed to link GGchem directly instead of shelling
  out, the calling convention it would need to match is unchanged by the
  version gap.

**Changed — the actual chemistry solver (`smchem16.f` / `smchem8.f`):**
This is where the two versions diverge the most. Current upstream adds:
- Explicit tracking of sulfur species as their own variables, mirroring the
  existing special-cased treatment of H2O/CO2/CH4:
  `pSO2,pSO3,pS2,pS8` and `lpSO2,lpSO3,lpS2,lpS8`, plus `SO2,SO3,SS2,SS8`
  molecule-index lookups via `stindex`.
- A dedicated correction/recovery path for COS (carbonyl sulfide):
  `COScorr`, `COSco`, `COSproblem`, alongside the pre-existing
  `HCOcorr`/`HCOproblem` mechanism.
- Numerical robustness guards not present in the vendored copy: a `lnpmin`
  clamp, a `sum>300*l` early-exit (`cycle`), and an `imaj2(i2)>0` guard before
  a particular branch is taken.
- A one-argument-signature change to the equilibrium-constant function call
  (`g(i)=gk(i)` → `g(i)=gk(i,Tg)`), and slightly more verbose diagnostic
  printing at `verbose>1`.

This amounts to a **solver robustness/convergence improvement specifically
for sulfur-bearing gas chemistry**, not the introduction of new molecules —
the reaction data itself (`dispol_BarklemCollet.dat`,
`dispol_StockKitzmann_withoutTsuji.dat`, `dispol_WoitkeRefit.dat`, as
actually read at runtime from `MARCS_CPH/MARCS/data/`) already lists
SO2/SO3/S2/S8 reactions in **both** the vendored and current data files
(counts match exactly except `dispol_StockKitzmann_withoutTsuji.dat`, where
upstream gained 4 extra lines — likely one additional reaction fit). Given
the ongoing Venus/Mars sulfur-chemistry gap work (OCS, Cl/HCl, S2, HF), this
is the one part of the version gap most likely to matter: the current
upstream solver should converge more reliably in S-rich regimes than the
vendored one, for the same input data.

**Changed — configuration surface (`read_parameter.f` / `datamod.f`):**
Upstream added several new optional run-time switches: `use_SiO`,
`metal_sulphates`, `output_dispol`, `model_refine`, `model_smooth`,
`Tmin_atmos`, `disk_model`, `adapt_file`. The `.in`-file parser is
keyword-tag based (`index(line,"! keyword")>0`, with defaults assigned
before parsing begins), so this is **backward compatible** — a
`marcs2ggchem.in` that doesn't mention these keys (as MARCS's does not)
simply gets upstream's defaults (`metal_sulphates=.true.`,
`use_SiO=.true.`, others off) rather than a deliberate, visible choice. No
format breakage; just inert new features until someone opts in.

**Changed — `init_chemistry.f`:** a small fix to how molecule names are
upper-cased before lookup (`uname=molname; call upper(uname)` replacing an
inline call), plus an optional `dispol.out` debug dump of parsed
equilibrium-constant fit coefficients (gated behind `output_dispol`, off by
default). Low risk either way.

**Changed — `equil_cond.f` / `supersat.f` (dust/condensation):** this is
where the upstream git log shows the most activity (condensate on/off
switching heuristics, H2SO4-related tweaks). **Currently inert for MARCS**:
the input file marcs.f writes hardcodes `model_eqcond = .false.`
(`marcs2ggchem_template.in`), so GGchem's equilibrium-condensation path is
never invoked from marcs — these changes have no effect on current MARCS
results regardless of which GGchem version is linked.

**Other differences** (`main.f`, `database.f`, `demo_*.f`, `auto_structure.f`,
`adapt_condensates.f`, the makefiles, and several files present in one tree
but not the other such as `disk_problem.f`, `calc_opac.f`, `MIEX.f90`,
`kappaRoss.f90`) belong to GGchem's standalone-driver/plotting/opacity
machinery, not the gas-phase chemistry path marcs.f exercises through
`call system(...)` — not relevant to marcs' results.

## Practical impact on MARCS today

- No input/output format breakage from upgrading: the `.in` parser is
  tolerant of missing keys, and `GGCHEM`'s call contract is unchanged.
- The one change that would alter actual chemistry output is the sulfur/COS
  solver robustness work in `smchem16.f`/`smchem8.f` — expect it to affect
  convergence/values specifically for S-bearing partial pressures, which
  matters for the Venus/Mars work.
- The condensation-related upstream changes are moot under the current
  MARCS call (`model_eqcond=.false.` always).
- The vendored `ggchem16` binary and source are functionally fine for
  non-sulfur chemistry; sulfur-focused work is the concrete reason to
  consider updating the vendored copy.

## Recommendation

Updating the vendored `MARCS_CPH/MARCS/GGchem/src16` copy to current
upstream looks safe from a plumbing standpoint (same call contract, tolerant
parser), but it **will** change chemistry results for sulfur-bearing layers
because of the smchem16/smchem8 solver changes. That means it should be
validated against a known/reference MARCS model (compare converged partial
pressures for S-bearing species before/after) rather than swapped in
silently. This is a separate follow-up decision, not part of the current
change (which only touched how `marcs.f` writes GGchem's input file, not
which GGchem version is linked).
