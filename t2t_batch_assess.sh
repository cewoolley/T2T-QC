#!/usr/bin/env bash
set -euo pipefail

# Reproduction of Monika Čechová's WDL pipeline using mamba envs
# Runs across all "*_assembly.bp.hap{1,2}.p_ctg.fa" assemblies and aggregates counts.
#
# Envs expected:
#   - t2t-assess : seqtk, bioawk, mashmap, gawk, sed, awk, paste, coreutils
#   - ncrf-py2  : NCRF, ncrf_summary
#
# Usage:
#   ./t2t_batch_assess_wdl.sh [-r ROOT] [-R REF] [-t THREADS]
#   ROOT defaults to ./spike_in_results next to this script
#   REF  defaults to ./hs1.fa next to this script
#   THREADS defaults to 8 (WDL default was 1; we allow override)
#
# Outputs (per assembly, identical names/columns as original wdl):
#   <name>.formatted.fa
#   <name>.edges.txt
#   <name>.headers.txt
#   <name>.unknown.txt
#   <name>.lengths.txt
#   <name>.combined.fasta
#   <name>.telomeric.summary.txt
#   <name>.telomeric.ends.txt
#   <name>.mashmap.txt
#   <name>.SUMMARY.txt
#   <name>.T2T.scaffolds.txt
#   <name>.T2T.contigs.txt
#
# + Aggregates outputs togehter into a final summary for plotting etc: 
#   ${ROOT}/t2t_summary.tsv   (variation  haplotype  t2t_scaffolds  t2t_contigs)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${ROOT:-${SCRIPT_DIR}/spike_in_results}"
REF_DEFAULT="${SCRIPT_DIR}/hs1.fa"
REF="${REF:-$REF_DEFAULT}"
THREADS="${THREADS:-8}"
BASES=1000         # WDL uses 1000 bases
TELO="TTAGGG"      # WDL searches TTAGGG only
AGG="${ROOT}/t2t_summary.tsv"

# ---- NEW: desired variation order (only used to prioritize) ----
ORDER=(pure_lsk ulk_01pct ulk_05pct ulk_10pct ulk_15pct ulk_20pct ulk_25pct ulk_30pct ulk_50pct pure_ulk)

usage() {
  echo "Usage: $0 [-r ROOT] [-R REF] [-t THREADS]"
  exit 1
}

while getopts ":r:R:t:" opt; do
  case $opt in
    r) ROOT="$OPTARG" ;;
    R) REF="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    *) usage ;;
  esac
done

[[ -d "$ROOT" ]] || { echo "ERROR: ROOT not found: $ROOT"; exit 1; }
[[ -f "$REF"  ]] || { echo "ERROR: REF FASTA not found: $REF"; exit 1; }

need() { command -v "$1" &>/dev/null || { echo "ERROR: required tool not found: $1"; exit 1; }; }
for t in mamba awk gawk sed paste seqtk bioawk mashmap; do need "$t"; done

# Activate mamba hooks
eval "$(mamba shell hook --shell bash)"

# --- helpers -----------------------------------------------------

# return 0 if $1 is contained in the rest of the args
contains() {
  local needle="$1"; shift || true
  local e
  for e in "$@"; do [[ "$e" == "$needle" ]] && return 0; done
  return 1
}

# Compute and append an AGG line for a given hap FASTA (must already have outputs)
write_agg_line() {
  local ASM="$1"
  local DIR NAME VAR HAP
  DIR="$(cd "$(dirname "$ASM")" && pwd)"
  NAME="$(basename "$ASM" .fa)"
  VAR="$(basename "$DIR")"
  if [[ "$NAME" =~ hap([12]) ]]; then HAP="hap${BASH_REMATCH[1]}"; else HAP="unknown"; fi
  local contigs_count scaffolds_count
  contigs_count=$( (tail -n +2 "${DIR}/${NAME}.T2T.contigs.txt" 2>/dev/null | awk 'NF>0' | wc -l) || echo 0 )
  scaffolds_count=$( (tail -n +2 "${DIR}/${NAME}.T2T.scaffolds.txt" 2>/dev/null | awk 'NF>0' | wc -l) || echo 0 )
  printf "%s\t%s\t%s\t%s\n" "$VAR" "$HAP" "$scaffolds_count" "$contigs_count" >> "$AGG"
}

# --- core worker (WDL-faithful) ---------------------------------
assess_one_faithful() {
  local ASM="$1"        # path to hap1/hap2 p_ctg.fa
  local REF="$2"        # reference fasta

  local DIR NAME VAR HAP
  DIR="$(cd "$(dirname "$ASM")" && pwd)"
  NAME="$(basename "$ASM" .fa)"                 # e.g., ulk_01pct_assembly.bp.hap1.p_ctg
  VAR="$(basename "$DIR")"                      # variation = parent directory
  if [[ "$NAME" =~ hap([12]) ]]; then HAP="hap${BASH_REMATCH[1]}"; else HAP="unknown"; fi

  pushd "$DIR" >/dev/null

  echo "### Running WDL-faithful assessment → ${VAR} | ${HAP} | ${NAME}"

  # --- formatAssembly (seqtk) ---
  mamba activate t2t-assess
  seqtk seq -l0 "$ASM" > "${NAME}.formatted.fa"

  # --- generateAssemblyEdges (gawk) ---
  # WDL assumes single-line fasta (header/sequence alternating lines).
  cat "${NAME}.formatted.fa" \
    | gawk -v len="$BASES" -F "" '{
        if (NR % 2 == 0) {
          for (i=1; i<=NF; i++) {
            if (i <= len || (NF-i) < len) { printf $(i) }
          }
          printf "\n"
        }
      }' > "${NAME}.edges.tmp.txt"
  cat "${NAME}.edges.tmp.txt" \
    | sed -e "s/.\{$BASES\}/&\n/g" | sed '/^$/d' > "${NAME}.edges.txt"
  rm -f "${NAME}.edges.tmp.txt"

  # --- runBioawk (unknown Ns, headers, lengths) ---
  bioawk -c fastx '{n=gsub(/N/, "", $seq); print $name "\t" n}' "$ASM" > "${NAME}.unknown.txt"
  bioawk -c fastx '{print ">"$name"_start";print ">"$name"_end";}' "$ASM" > "${NAME}.headers.txt"
  bioawk -c fastx '{print $name, length($seq)}' "$ASM" > "${NAME}.lengths.txt"

  # --- combineIntoFasta (paste) ---
  paste -d $'\n' "${NAME}.headers.txt" "${NAME}.edges.txt" > "${NAME}.combined.fasta"

  # --- runNCRF (TTAGGG only, matches WDL flags) ---
  mamba deactivate
  mamba activate ncrf-py2
  need NCRF
  need ncrf_summary
  cat "${NAME}.combined.fasta" \
    | NCRF --stats=events "$TELO" \
    | ncrf_summary > "${NAME}.telomeric.summary.txt"
  tail -n +2 "${NAME}.telomeric.summary.txt" | cut -f3 | sort | uniq > "${NAME}.telomeric.ends.txt"

  # --- runMashMap ---
  mamba deactivate
  mamba activate t2t-assess
  mashmap --threads "$THREADS" --perc_identity 95 --noSplit \
    -r "$REF" -q "$ASM" -o "${NAME}.mashmap.txt"

  # --- assessCompleteness (robust to missing matches; safe under -euo pipefail) ---
  echo "Task 7: assessCompleteness (robust)"
  {
  while IFS= read -r line; do
    contig="${line%%_*}"  # strip _start/_end

    # lengths: contig<space>length (fallback to contig<TAB>NA)
    { grep -m1 -F -- "$contig" "${NAME}.lengths.txt" || printf '%s\tNA\n' "$contig"; } | tr -d '\n'
    printf '\t'

    # unknown Ns: contig<TAB>Ns (fallback to contig<TAB>NA)
    { grep -m1 -F -- "$contig" "${NAME}.unknown.txt" || printf '%s\tNA\n' "$contig"; } | tr -d '\n'
    printf '\t'

    # mashmap PAF line for contig (first hit). If none, emit a stub with >=12 fields so cut -f10 exists.
    { grep -m1 -F -- "$contig" "${NAME}.mashmap.txt" \
      || printf '%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n' "$contig"; } \
      | tr -s '[:blank:]' '\t'

  done < "${NAME}.telomeric.ends.txt"
  } | cut -f1,2,4,10 \
    | sort | uniq -c | sort -rgk1 \
    | tr -s '[:blank:]' '\t' | sed 's/^\t//' \
    > "${NAME}.SUMMARY.txt"

  # Output tables with headers exactly as WDL
  echo -e 'name\tlength\tNs\tchromosome' > "${NAME}.T2T.scaffolds.txt"
  echo -e 'name\tlength\tNs\tchromosome' > "${NAME}.T2T.contigs.txt"

  if grep -q '^2' "${NAME}.SUMMARY.txt"; then
    grep '^2' "${NAME}.SUMMARY.txt" | awk '{if ($4>=0) print;}' \
      | cut -f2- | sort -k4,4 -V -s >> "${NAME}.T2T.scaffolds.txt" || true
    grep '^2' "${NAME}.SUMMARY.txt" | awk '{if ($4==0) print;}' \
      | cut -f2- | sort -k4,4 -V -s >> "${NAME}.T2T.contigs.txt" || true
  else
    echo "None of the contigs/scaffolds contained telomeres on both ends."
  fi

  # Print a quick per-hap summary to stdout (aggregator handled outside)
  local contigs_count scaffolds_count
  contigs_count=$( (tail -n +2 "${NAME}.T2T.contigs.txt" | awk 'NF>0' | wc -l) || echo 0 )
  scaffolds_count=$( (tail -n +2 "${NAME}.T2T.scaffolds.txt" | awk 'NF>0' | wc -l) || echo 0 )

  echo "✓ ${VAR} | ${HAP} | scaffolds=${scaffolds_count} contigs=${contigs_count}"
  popd >/dev/null
}

# --- discovery ---------------------------------------------------
echo "Scanning for *_assembly.bp.hap{1,2}.p_ctg.fa under: $ROOT"
# shellcheck disable=SC2207
readarray -d '' FASTAS < <(find "$ROOT" -type f -name "*_assembly.bp.hap[12].p_ctg.fa" -print0)

if (( ${#FASTAS[@]} == 0 )); then
  echo "No hap1/hap2 p_ctg FASTA files found. Run the GFA→FA conversion first."
  exit 1
fi

# Map (variation,hap) -> FASTA; collect variations seen
declare -A FASTA_BY_VAR_HAP
declare -A VARS_SEEN
for fa in "${FASTAS[@]}"; do
  DIR="$(cd "$(dirname "$fa")" && pwd)"
  VAR="$(basename "$DIR")"
  BASE="$(basename "$fa" .fa)"
  HAP="unknown"
  if [[ "$BASE" =~ hap([12]) ]]; then HAP="hap${BASH_REMATCH[1]}"; fi
  FASTA_BY_VAR_HAP["$VAR|$HAP"]="$fa"
  VARS_SEEN["$VAR"]=1
done

# Print the variations discovered (unique)
echo
echo "Found variations (${#VARS_SEEN[@]}):"
printf '  - %s\n' $(printf '%s\n' "${!VARS_SEEN[@]}" | sort)

# Build processing order: requested ORDER first (if present), then extras sorted
ordered_present=()
for v in "${ORDER[@]}"; do
  if [[ -n "${VARS_SEEN[$v]+x}" ]]; then ordered_present+=( "$v" ); fi
done
extras=()
for v in "${!VARS_SEEN[@]}"; do
  contains "$v" "${ordered_present[@]}" || extras+=( "$v" )
done
# sort extras lexicographically
IFS=$'\n' extras_sorted=($(printf "%s\n" "${extras[@]:-}" | sort || true)); unset IFS
final_vars=( "${ordered_present[@]}" "${extras_sorted[@]:-}" )

echo
echo "Processing order:"
printf '  %2d. %s\n' $(i=1; for v in "${final_vars[@]}"; do echo "$((i++))" "$v"; done)

# --- aggregator header (rebuilt fresh every run) -----------------
printf "variation\thaplotype\tt2t_scaffolds\tt2t_contigs\n" > "$AGG"

# --- execution with resume flags --------------------------------
for VAR in "${final_vars[@]}"; do
  for HAP in hap1 hap2; do
    ASM="${FASTA_BY_VAR_HAP["$VAR|$HAP"]:-}"
    if [[ -z "${ASM}" ]]; then
      echo "[$VAR][$HAP] FASTA not found — skipping."
      continue
    fi

    DIR="$(cd "$(dirname "$ASM")" && pwd)"
    BASE="$(basename "$ASM" .fa)"
    DONE_FLAG="${DIR}/.${BASE}.done"
    RUN_FLAG="${DIR}/.${BASE}.running"

    # If outputs exist but no flag (e.g., first run with this resume logic), mark done and skip
    if [[ ! -f "$DONE_FLAG" ]] && [[ -f "${DIR}/${BASE}.T2T.contigs.txt" ]] && [[ -f "${DIR}/${BASE}.T2T.scaffolds.txt" ]]; then
      echo "[$VAR][$HAP] Outputs present but no .done flag — marking as done."
      touch "$DONE_FLAG"
    fi

    if [[ -f "$DONE_FLAG" ]]; then
      echo "[$VAR][$HAP] already completed — skipping."
      write_agg_line "$ASM"
      continue
    fi

    # Clean any stale .running flag
    [[ -f "$RUN_FLAG" ]] && rm -f "$RUN_FLAG"
    touch "$RUN_FLAG"

    # Run the assessment
    assess_one_faithful "$ASM" "$REF"

    # Mark success and update aggregator; then clear running flag
    rm -f "$RUN_FLAG"
    touch "$DONE_FLAG"
    write_agg_line "$ASM"
  done

done

echo
echo "Summary written to: $AGG"
# Pretty print, fallback to plain if 'column' not available
if command -v column >/dev/null 2>&1; then
  column -t -s $'\t' "$AGG" || cat "$AGG"
else
  cat "$AGG"
fi
