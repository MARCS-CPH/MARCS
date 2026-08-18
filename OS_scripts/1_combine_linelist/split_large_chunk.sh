#!/bin/bash
# Split an oversized ExoMol .trans chunk into ~24G line-safe pieces and submit
# one long-running job that processes each piece sequentially through
# loop_v2.py, skipping any piece that already has a line1.txt.
#
# When you need this: a single 1000cm-1 (or however wide) linelist chunk is
# too big to run through loop_v2.py directly -- either it OOMs at moderate
# memory, or a plain run silently "completes" with no line1.txt because
# multiprocessing swallows a crashed child (loop_v2.py has no exitcode
# check -- always verify line1.txt actually exists, don't trust SLURM
# COMPLETED status alone). Splitting into smaller pieces and processing them
# one at a time in a single long job is the proven workaround.
#
# Usage: split_large_chunk.sh <SPECIES> <EXOMOL_TAG> <CHUNK> <TIME_LIMIT>
#   SPECIES     e.g. H2CS        (matches /lustre/astro/tbalduin/molecules/<SPECIES>/)
#   EXOMOL_TAG  e.g. 1H2-12C-32S__MOTY   (the raw ExoMol linelist basename --
#               check linelist/*.states to find it for your species)
#   CHUNK       e.g. 06000-07000   (the wavenumber-range chunk to split)
#   TIME_LIMIT  e.g. 14:00:00      (generous -- ~30min/24G piece is typical,
#               so pieces * 30min + margin; use astro3_long, not
#               astro3_short, since that partition caps at 12h)
#
# Example: 359GB chunk -> ~14 pieces * ~30min = ~7h -> use 12:00:00
set -euo pipefail

SPECIES="$1"
TAG="$2"
CHUNK="$3"
TIME_LIMIT="$4"

BASE="/lustre/astro/tbalduin/molecules/${SPECIES}/linelist"
PARENT="${BASE}/${TAG}__${CHUNK}"
TRANS="${PARENT}/${TAG}__${CHUNK}.trans"
STATES_LINK="${PARENT}/${TAG}.states"

echo "=== Splitting $TRANS into ~24G line-safe pieces ==="
split -C 24G -d --numeric-suffixes=1 -a 3 --additional-suffix=.trans_raw "$TRANS" "${PARENT}/splitpart_"

n=0
for f in "${PARENT}"/splitpart_*.trans_raw; do
  n=$((n+1))
  idx=$(printf "%03d" "$n")
  d="${BASE}/${TAG}__${CHUNK}_splitpart_${idx}"
  mkdir -p "$d"
  mv "$f" "${d}/splitpart_${idx}.trans"
  ln -sf "$STATES_LINK" "${d}/${TAG}.states"
  cp "${PARENT}/loop_v2.py" "${d}/loop_v2.py"
done
echo "Created $n splitpart directories."

DRIVER="${PARENT}/run_splits.sh"
cat > "$DRIVER" << EOF
#!/bin/bash
#SBATCH -p astro3_long
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH -N 1
#SBATCH --mem=650gb
#SBATCH --time=${TIME_LIMIT}
#SBATCH --job-name=${SPECIES,,}_split_${CHUNK}
#SBATCH --output=${PARENT}/slurm-%j.out

for d in ${BASE}/${TAG}__${CHUNK}_splitpart_*/; do
  echo "=== \$(date) starting \$d ==="
  cd "\$d"
  if [ -f line1.txt ]; then
    echo "already has line1.txt, skipping"
    continue
  fi
  python loop_v2.py >> loop.runinfo.out
  if [ -f line1.txt ]; then
    echo "=== \$(date) OK: \$d ==="
  else
    echo "=== \$(date) FAILED (no line1.txt produced): \$d ==="
  fi
done
echo "=== \$(date) ALL SPLITPARTS DONE for ${CHUNK} ==="
EOF

echo "=== Submitting driver job (time limit ${TIME_LIMIT}) ==="
sbatch "$DRIVER"
