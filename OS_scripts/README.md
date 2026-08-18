# OS_scripts — MARCS opacity-sampling (OS) file build pipeline

Scripts used to build new molecular opacity-sampling files for MARCS
(`data/mol_names_cool_75.dat` for example), covering both the full-linelist path
(ExoMol download → `os26`) and the UV-cross-section-only fallback path for
species with no rovibrational linelist available anywhere. Built out and
debugged over the course of adding ~45 new species to this project's WIP
opacity list; every gotcha noted below actually happened at least once.

## Paths used below

The commands and script defaults reference a few external directories that
are site-specific — swap in your own locations:

| Placeholder | What it is |
|---|---|
| `<PATH TO ORIGINAL LINELIST>` | Root directory holding one subdirectory per species' raw ExoMol download + `os26` working area (linelist chunks, `.pf` file, `crosssections/`) |
| `<PATH TO OS26 TOOL>` | Directory containing the compiled `os26` binary + source |
| `<PATH TO OPACITY CONVERTER>` | The `marcs_opac_converter` package checkout (stage 3, `convert_opac`) |
| `<PATH TO UV CROSSSECTIONS>` | Root directory holding literature UV cross-section archives (Leiden/SWRI/MPI-Mainz format subfolders) |

The scripts in this folder hardcode this project's actual paths as
defaults (e.g. `CROSSSECS = '/groups/.../crosssecs'` at the top of the
Python scripts, `BASE="/lustre/.../molecules"` in the shell scripts) —
edit those constants for your own layout before reusing them elsewhere.

## Two build paths

**Path A — full rovibrational linelist exists (ExoMol/HITRAN).** Produces
real temperature- and wavenumber-resolved opacity. Use `1_combine_linelist/`
→ `2_os26_crosssections/` → `3_convert_opac/`.

**Path B — no rovibrational linelist, only a literature UV cross section
(Leiden/SWRI/MPI-Mainz).** UV-only opacity, no IR branch. Use
`4_uvonly_alternative/` → `3_convert_opac/`. Much faster (minutes, not
days) since there's no giant linelist to process — the whole cost is stage
3's `os26` run in Path A.

Both paths converge on stage 3 (`convert_opac`) and stage 5 (wiring into
`mol_names`).

```
1_combine_linelist/     Path A, stage 1-2: combine chunked ExoMol linelist
                         downloads into one master line1.txt per species
2_os26_crosssections/   Path A, stage 3: run os26 to sample cross sections
                         across the MARCS temperature/wavenumber grid
3_convert_opac/         Both paths, stage 4: convert to the binary
                         crossec.dat/wn.dat/input.nml triplet MARCS reads
4_uvonly_alternative/   Path B: build a UV-only intermediate file directly
                         from a literature cross section, skipping 1 and 2
5_wire_species/         Both paths, stage 5: add the species to a
                         mol_names file
```

## Quick start

**Path A (full linelist), once the raw ExoMol chunks are downloaded to
`<PATH TO ORIGINAL LINELIST>/<SPECIES>/linelist/`:**
```bash
# 1. If any single chunk is too big to run through loop_v2.py directly
#    (OOMs, or silently produces no output despite SLURM reporting
#    COMPLETED -- always check, see script comments):
1_combine_linelist/split_large_chunk.sh <SPECIES> <EXOMOL_TAG> <CHUNK> <TIME_LIMIT>

# 2. Verify every chunk has a real line1.txt before combining
1_combine_linelist/verify_chunks_complete.sh <SPECIES> <EXOMOL_TAG>

# 3. Combine into one master line1.txt
1_combine_linelist/combine_chunks.sh <SPECIES> <EXOMOL_TAG>

# 4. Set up and run os26 (see 2_os26_crosssections/ template comments,
#    and worked_example_C2H4/ for a complete real example)
mkdir -p <PATH TO ORIGINAL LINELIST>/<SPECIES>/crosssections
ln -s <PATH TO OS26 TOOL>/os26 \
      <PATH TO ORIGINAL LINELIST>/<SPECIES>/crosssections/os26
ln -s <PATH TO ORIGINAL LINELIST>/<SPECIES>/linelist/<TAG>.pf \
      <PATH TO ORIGINAL LINELIST>/<SPECIES>/linelist/<MOLID>.pf
# fill in os.input from the template, checking TMOL against the .pf file's
# real range (see template comments -- this is the #1 way this breaks)
cp 2_os26_crosssections/os.input.template \
   <PATH TO ORIGINAL LINELIST>/<SPECIES>/crosssections/os.input
cp 2_os26_crosssections/run_full.sh.template \
   <PATH TO ORIGINAL LINELIST>/<SPECIES>/crosssections/run_full.sh
cd <PATH TO ORIGINAL LINELIST>/<SPECIES>/crosssections && sbatch run_full.sh

# 5. Once done, verify before trusting it
grep -c Infinity full_run/<SPECIES>_OS_*.dat    # must print 0
```

**Path B (UV-only), once you've confirmed the species has no linelist but
does have a local/MPI-Mainz UV cross section:**
```bash
# Edit the SPECIES list in build_uvonly_flat.py (or build_uvonly_multiT.py
# if a genuine multi-temperature series exists), then:
python3 4_uvonly_alternative/build_uvonly_flat.py
```

**Both paths converge here:**
```bash
# stage 3: convert to the MARCS binary triplet
cp 3_convert_opac/config_stage4.yaml.template config_stage4_<species>.yaml
# edit it, then:
cd <PATH TO OPACITY CONVERTER>
python3 convert_opac config_stage4_<species>.yaml

# stage 5: wire into a mol_names file (WIP copy, not the live one, unless
# you're sure nothing is currently reading it)
5_wire_species/add_species.sh <MOL_NAMES_FILE> <KEY> <OUTPUT_DIR>
```

## Gotchas (all encountered for real, not hypothetical)

- **`loop_v2.py` (linelist chunk processing) can silently fail.** It uses
  `multiprocessing.Process().join()` with no exit-code check, so a crashed
  or OOM-killed child doesn't fail the parent job — SLURM reports
  `COMPLETED`, exit code `0`, and there's simply no `line1.txt` written.
  Also seen: a job that times out mid-run leaves nothing behind either,
  same silent-looking failure. **Always verify `line1.txt` exists**
  (`verify_chunks_complete.sh`), never trust job state alone.

- **`os26`'s TMOL grid must not exceed the partition-function file's own
  temperature range**, or every column beyond that range comes out as the
  literal string `Infinity` — no warning, no error, it just writes it.
  Different species' ExoMol-downloaded `.pf` files cover different native
  ranges (not related to the molecule's real thermal stability — purely an
  artifact of what that particular ExoMol database entry ships). Check
  `tail -3 linelist/<TAG>.pf` before picking TMOL. Don't infer a T-range
  "by analogy" to a chemically similar species — check the actual file.

- **`os26` needs `linelist/<MOLID>.pf`, not `linelist/<TAG>.pf`.** The raw
  ExoMol download names the partition-function file after its own dataset
  tag (e.g. `12C2-1H4__MaYTY.pf`), but `os26` looks for it under the
  MOLID you put in `os.input` (e.g. `C2H4.pf`). Missing symlink = instant
  crash (`forrtl: severe (29): file not found, unit 67`), 0-second job,
  easy to miss if you're not watching closely.

- **SWRI cross-section files are not simple fixed-column tables.** They
  concatenate multiple reference/branching-ratio blocks in one file, each
  with its own free-text header of varying length, followed by a numeric
  table. Read only the FIRST contiguous numeric block (the "Total" or sole
  first channel) and stop — don't blindly skip N header lines and read to
  EOF, you'll pull in a second, differently-shaped table. Also: the total
  cross-section column position varies per file, always check it directly
  rather than reusing another species' hardcoded column index.

- **Wavenumber-to-wavelength unit conversion reverses sort order.**
  `wavenumber = 1e8 / wavelength_Å` is order-reversing (1/x). Sort AFTER
  converting, never before — sorting wavelength ascending then converting
  silently produces a fully-reversed (descending) wavenumber array, which
  MARCS refuses to load. This actually happened across all three Path-B
  build scripts in one pass (copy-pasted the same bug three times) before
  being caught by a downstream MARCS test-run bisection — see the fix
  history in `build_uvonly_flat.py`'s docstring.

- **Runtime for `os26` scales with the combined `line1.txt` size, roughly
  linearly** (~0.029-0.030 h/GB observed on this cluster). Don't leave
  `--time` at a short default for anything over a few hundred GB — use
  a long-running partition (multi-day cap) and size generously. A long
  final sort/write phase after the input-reading phase finishes can add
  many more hours on top of when the scratch file (`lines.tmp`) stops
  growing — that's normal, not a hang, for the largest species.

- **Lustre `du` can under-report by ~2x for very recently written large
  files** — its block-accounting has real propagation lag (matches the
  same lag `lfs quota` shows right after a big write or delete). Don't
  trust a `du -sh` reading taken immediately after a large operation
  completes; wait a few seconds and recheck, or use `stat -c%s` /
  `wc -c` for an authoritative byte count.

- **Disk quota has a grace period, not a hard cutoff** — exceeding the
  soft quota starts roughly a week-long grace timer before writes actually
  start failing, it doesn't block immediately. Don't panic-clean under
  time pressure; check `lfs quota -u <user> <filesystem>` for the actual
  grace countdown before treating it as an emergency.

- **What's actually safe to delete, in order of "how much it costs to
  regenerate" (cheapest to most expensive):**
  1. `os26`'s `lines.tmp` scratch file (`crosssections/full_run/`) — pure
     scratch, safe to delete the moment the corresponding `<SPECIES>_OS_*.dat`
     exists and is verified (no `Infinity`, spot-checked a data row).
  2. Raw `.trans` ExoMol download chunks, once every chunk's `line1.txt`
     is verified present — these are the original downloaded data, cheap
     to re-fetch if ever needed again.
  3. Per-chunk `line1.txt` files, once the combined master `line1.txt` is
     verified byte-identical to their sum — these are the expensive part
     to redo (real compute, hours per large chunk, and the source `.trans`
     data would already be gone by this point in the cleanup order), so
     prefer freeing `.trans` first if space alone is the goal.

- **A species' left-column key in `mol_names` doesn't have to match its
  OS file's internal `MOLID`.** Established convention in this project:
  `HS` (network key) → `SH_new/` (dir, internal `MOLID='SH'`); `HNO2` →
  `HONO_new_flatT/` (internal `MOLID='HONO'`). This is intentional, not a
  bug to reconcile — the network/KROME token and the spectroscopic
  database's own naming can legitimately differ.

- **Never edit the live `mol_names_cool_75.dat` while a MARCS run might be
  reading it.** Stage new species in a separate WIP copy
  (`mol_names_cool_WIP.dat` in this project) and only merge into the live
  file when you're sure nothing is using it, or ask first.

## Reference locations (not part of this repo, linked for context)

- Raw ExoMol linelist downloads + `os26` working area:
  `<PATH TO ORIGINAL LINELIST>/<SPECIES>/`
- `os26` binary + source: `<PATH TO OS26 TOOL>/`
- Literature UV cross sections (Leiden/SWRI/MPI-Mainz):
  `<PATH TO UV CROSSSECTIONS>/{Leiden,SWRI,MPI}/`
- `marcs_opac_converter` package (stage 3, `convert_opac`):
  `<PATH TO OPACITY CONVERTER>` — not pip-installed, only importable when
  run from inside that directory (relies on Python's implicit CWD
  sys.path entry)
- Species-gap / temperature-coverage tracking docs, if you keep similar
  ones for your own project: useful to record which species still need an
  OS file and why, and which existing ones are on a stale temperature grid
