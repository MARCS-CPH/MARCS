#!/bin/bash
# Check that every linelist chunk for a species has a genuine line1.txt --
# either directly, or via all of its splitpart_NNN/ subdirectories if it was
# split (see split_large_chunk.sh). Run this before combine_chunks.sh, and
# again before deleting any .trans files.
#
# Why this matters: loop_v2.py uses multiprocessing.Process().join() with no
# exitcode check, so a child that OOMs or crashes doesn't fail the parent --
# the job can report SLURM COMPLETED / exit 0 while silently producing no
# line1.txt at all for one or more chunks. Also seen: a job that times out
# mid-run leaves a partial/absent line1.txt with no error either. Always
# verify the actual file, never trust job exit status alone.
#
# Usage: verify_chunks_complete.sh <SPECIES> <EXOMOL_TAG>
set -euo pipefail

SPECIES="$1"
TAG="$2"
BASE="/lustre/astro/tbalduin/molecules/${SPECIES}/linelist"

missing=0
checked=0
for d in $(ls -d "${BASE}/${TAG}__"*/ 2>/dev/null | grep -v splitpart | sort); do
  checked=$((checked + 1))
  chunk=$(basename "$d" | sed "s/${TAG}__//")
  if [ -f "${d}line1.txt" ]; then
    continue
  fi
  splitparts=$(ls -d "${d%/}_splitpart_"*/ 2>/dev/null || true)
  if [ -z "$splitparts" ]; then
    echo "MISSING: $chunk has no line1.txt and no splitpart subdirs"
    missing=$((missing + 1))
    continue
  fi
  for sp in $splitparts; do
    if [ ! -f "${sp}line1.txt" ]; then
      echo "MISSING: $chunk splitpart $(basename "$sp") has no line1.txt"
      missing=$((missing + 1))
    fi
  done
done

echo "=== checked $checked chunks, $missing incomplete ==="
[ "$missing" -eq 0 ]
