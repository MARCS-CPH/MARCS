#!/bin/bash
# Combine all of a species' linelist chunks into one master line1.txt,
# auto-discovering chunk directories under linelist/. Handles both plain
# chunks (a direct line1.txt in the chunk dir) and chunks that were split
# via split_large_chunk.sh (multiple <chunk>_splitpart_NNN/ subdirs, each
# with its own line1.txt -- concatenated in numeric order).
#
# Usage: combine_chunks.sh <SPECIES> <EXOMOL_TAG>
#   SPECIES     e.g. C2H4, H2CS
#   EXOMOL_TAG  e.g. 12C2-1H4__MaYTY   (basename before the __<chunk> suffix
#               -- check linelist/*.states to find it)
#
# Before running: every chunk directory must already have a verified
# line1.txt (direct, or via all its splitpart_NNN/ subdirs) -- os26 will
# happily read a partial/incomplete line1.txt with no error, so check
# completeness first rather than trusting that the combine "worked" just
# because it produced a file. This script itself doesn't verify each
# source -- run 1_combine_linelist/verify_chunks_complete.sh first, or eyeball
# it: `for d in linelist/*/; do [ -f "$d/line1.txt" ] || echo "$d"; done`
# should print nothing once splitpart dirs are excluded from that check.
#
# Disk note: this creates a NEW file the same size as the sum of all chunk
# line1.txt files, without deleting the sources. For big species (multi-TB)
# check free space first -- once the combined master is verified (e.g. byte
# count matches: `du -shc linelist/*/line1.txt` sums to the master's size),
# the per-chunk line1.txt files and the raw .trans files are both safe to
# delete (the .trans files are pure ExoMol download data, cheap to
# re-fetch if ever needed again; line1.txt is the more expensive-to-redo
# asset -- prefer freeing .trans first if you need the space).
set -euo pipefail

SPECIES="$1"
TAG="$2"

BASE="/lustre/astro/tbalduin/molecules/${SPECIES}/linelist"
OUT="${BASE}/line1.txt.tmp"
> "$OUT"

# discover chunk dirs (exclude splitpart subdirs) in sorted wavenumber order
n_chunks=0
for d in $(ls -d "${BASE}/${TAG}__"*/ 2>/dev/null | grep -v splitpart | sort); do
  n_chunks=$((n_chunks + 1))
  chunk=$(basename "$d" | sed "s/${TAG}__//")
  if [ -f "${d}line1.txt" ]; then
    cat "${d}line1.txt" >> "$OUT"
    echo "chunk $chunk: direct line1.txt appended"
  else
    n=0
    for sp in $(ls -d "${d%/}_splitpart_"*/ 2>/dev/null | sort); do
      cat "${sp}line1.txt" >> "$OUT"
      n=$((n + 1))
    done
    if [ "$n" -eq 0 ]; then
      echo "chunk $chunk: NO line1.txt found (direct or split) -- ABORTING"
      rm -f "$OUT"
      exit 1
    fi
    echo "chunk $chunk: $n splitparts appended"
  fi
done

if [ "$n_chunks" -eq 0 ]; then
  echo "No chunk directories found matching ${BASE}/${TAG}__*/ -- check SPECIES/EXOMOL_TAG"
  rm -f "$OUT"
  exit 1
fi

mv "$OUT" "${BASE}/line1.txt"
echo "=== combined ${n_chunks} chunks ==="
wc -l "${BASE}/line1.txt"
stat -c '%s bytes' "${BASE}/line1.txt"
