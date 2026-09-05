# KROME chemical network versions

**Note on filenames vs. the names used here:** "v6" below is `krome/networks/non_eq/full_network_DMS_v2.dat` and "v7" is `full_network_v8.dat` — the on-disk names are historical and pre-date this simplification. The repo also still has an older, unrelated `full_network_v6.dat` (122 species/883 reactions) that is **not** the "v6" described here; it's superseded and kept only for compatibility with older scripts that reference it by that literal filename.

## The networks

- **v4** — earliest version, flat photo rate constants (no real cross-sections).
- **v5_full_photo** — adds real photolysis cross-sections (`@xsecFile=LEIDEN`).
- **v5_static_photo** — same base as v5, photochemistry stripped out.
- **v6** — adds the full organic/sulfur suite.
- **full_network_earth** — v6 with the CH2O/HCHO duplicate-species bug fixed. Same size class as v6, not v7.
- **v7** — adds more hydrocarbons on top of v6, plus a duplicate-species bugfix. The largest network.

## Size comparison

| Network | Species | Reactions (non-photo) | Reactions (photo) | Total reactions |
|---|---:|---:|---:|---:|
| v4 | 49 | 354 | 34 | 388 |
| v5_full_photo | 51 | 359 | 40 | 399 |
| v5_static_photo | 50 | 396 | 0 | 396 |
| v6 | 94 | 387 | 86 | 473 |
| full_network_earth | 93 | 387 | 86 | 473 |
| v7 | 137 | 843 | 98 | 941 |

Species counts for the built networks (v6, full_network_earth, v7) include 5 non-chemical bookkeeping placeholders (M, CR, g, Tgas, dummy); real chemical-species counts run ~5 lower.

## What each step adds

- **v4 → v5**: +2 species (`C`, `CH`).
- **v5 → v6** (large jump): the whole organic/sulfur suite arrives at once — C2-C6 hydrocarbons, the DMS/organosulfur branch (`CH3SCH3`=DMS, `CH3SH`, `CH3SO`, `CH3SOCH3`, `C5H10O2S`=DMSP, `C2H6SO`, `C4H10S(2)`, `COS`/`OCS`, `CS2`, `S3/S4/S8`), the nitrile/HCN family, and the N-H "hydrazine ladder" (`NH`/`NH2`/`NH3`, `N2H`/`N2H2`/`N2H3`/`N2H4`).
- **v6 → full_network_earth**: merges `CH2O` into `HCHO` — same molecule (formaldehyde), previously tracked as two disconnected species. Branches off v6 directly, not derived from v7.
- **v6 → v7** (separate branch from full_network_earth): adds more hydrocarbons (C2-C4 range) and fixes a naming-duplicate bug (`HN2`, a duplicate of `N2H`, removed).
