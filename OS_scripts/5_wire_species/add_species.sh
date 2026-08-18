#!/bin/bash
# Safely append one species entry to a mol_names file: checks the species
# isn't already present (by exact key match), checks the target directory
# has a complete crossec.dat/wn.dat/input.nml triplet, and handles the
# missing-trailing-newline gotcha (mol_names_cool_75.dat's last line has no
# trailing newline -- a naive `echo >> file` would glue your new entry onto
# the end of the previous line instead of starting a new one).
#
# Usage: add_species.sh <MOL_NAMES_FILE> <KEY> <OUTPUT_DIR>
#   MOL_NAMES_FILE  e.g. /lustre/hpc/astro/tbalduin/MARCS_CPH/MARCS/data/mol_names_cool_WIP.dat
#                   NEVER point this at the live mol_names_cool_75.dat while
#                   MARCS runs are using it -- stage new species in a WIP
#                   copy, merge into live only when you're ready to.
#   KEY             the left-column species identifier MARCS/KROME look up,
#                   e.g. C2H4. Doesn't have to match the directory's own
#                   internal MOLID -- see the HS/SH and HNO2/HONO examples
#                   in this project's mol_names_cool_WIP.dat: the network
#                   token and the OS file's own MOLID can legitimately
#                   differ, that's not a bug to "fix".
#   OUTPUT_DIR      the output/combined/<X>/ directory from stage 3, with a
#                   trailing slash (matches this project's existing
#                   convention in mol_names files).
#
# Example:
#   add_species.sh data/mol_names_cool_WIP.dat C2H4 \
#     /groups/astro/tbalduin/marcs_opac_converter/combined_UV_IR/output/combined/C2H4_new_100to1286K/
set -euo pipefail

MOL_NAMES="$1"
KEY="$2"
OUTDIR="$3"

if [ ! -f "$MOL_NAMES" ]; then
  echo "ERROR: $MOL_NAMES does not exist"
  exit 1
fi

if [ ! -f "${OUTDIR}crossec.dat" ] || [ ! -f "${OUTDIR}wn.dat" ] || [ ! -f "${OUTDIR}input.nml" ]; then
  echo "ERROR: $OUTDIR is missing crossec.dat/wn.dat/input.nml -- run stage 3 (convert_opac) first"
  exit 1
fi

if tail -n +2 "$MOL_NAMES" | awk '{print $1}' | grep -qx "$KEY"; then
  echo "ERROR: '$KEY' is already present in $MOL_NAMES -- not adding a duplicate"
  echo "  (existing line: $(grep "^${KEY}[[:space:]]" "$MOL_NAMES"))"
  exit 1
fi

# ensure the file ends with a newline before appending, or the new entry
# glues onto the last existing line
if [ -s "$MOL_NAMES" ] && [ "$(tail -c 1 "$MOL_NAMES")" != "" ]; then
  echo >> "$MOL_NAMES"
fi

echo "${KEY} ${OUTDIR}" >> "$MOL_NAMES"

echo "Added: ${KEY} -> ${OUTDIR}"
echo "New line count: $(wc -l < "$MOL_NAMES")"

# final duplicate-key sanity check across the whole file
dupes=$(tail -n +2 "$MOL_NAMES" | awk '{print $1}' | sort | uniq -d)
if [ -n "$dupes" ]; then
  echo "WARNING: duplicate keys now present in $MOL_NAMES:"
  echo "$dupes"
fi
