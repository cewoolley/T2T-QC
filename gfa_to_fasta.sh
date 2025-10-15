#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./gfa_to_fasta.sh [ROOT_OUTPUT_DIR]
#
# If not provided, ROOT defaults to "<this-script-dir>/spike_in_results"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${1:-${SCRIPT_DIR}/spike_in_results}"

echo "Scanning for hap1/hap2 p_ctg GFAs under: ${ROOT}"
echo

# Handle filenames safely (spaces, etc.)
converted=0
found_any=false

# Find hifiasm-ont partially phased contig GFAs and convert each to FASTA in-place
while IFS= read -r -d '' gfa; do
  found_any=true
  fa="${gfa%.gfa}.fa"
  echo "Converting:"
  echo "  GFA: ${gfa}"
  echo "   →  FA: ${fa}"

  # GFA is tab-delimited; fields: col2 = segment name, col3 = sequence
  # Using /^S\t/ guards we only emit 'S' records.
  awk 'BEGIN{FS="\t"} /^S\t/ {print ">" $2 "\n" $3}' "$gfa" > "$fa"
  converted=$((converted+1))
done < <(find "$ROOT" -type f -name "*_assembly.bp.hap[12].p_ctg.gfa" -print0)

if ! $found_any; then
  echo "No matching '*_assembly.bp.hap[12].p_ctg.gfa' files found under: ${ROOT}"
else
  echo
  echo "Done. Converted ${converted} GFA file(s) to FASTA."
fi
