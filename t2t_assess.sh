#!/bin/bash
set -euo pipefail

echo "### Mamba (Python) Re-implementation of Monika Cechova's T2T Assessment Workflow"
echo "### Connor Woolley"


# needs one arg
if [ "$#" -ne 1 ]; then
        echo "Usage: $0 <filename.fa>"
        exit 1
fi

eval "$(mamba shell hook --shell bash)"

echo "Moving to T2T conda env..."
mamba activate t2t-assess

#ASM="hap1.fa"      # Swap to hap2.fa for second run; now an arg, is better.
ASM=$1
REF="hs1.fa"
NAME=$(basename "$ASM" .fa)  # e.g., hap1

echo "Running T2T assessment for $ASM..."


echo "Task 1: formatAssembly (seqtk) - rewrap to single-line FASTA"
seqtk seq -l0 "$ASM" > "${NAME}.formatted.fa" || { echo "Task 1 failed"; exit 1; }

echo "Task 2: generateAssemblyEdges (gawk) - extract 1000 bp ends"
bases=1000
gawk -v len="$bases" -F '' 'NR % 2 == 0 {
  for (i=1; i<=NF; i++) {
    if (i <= len || (NF-i) < len) { printf $i }
  }
  printf "\n"
}' "${NAME}.formatted.fa" > "${NAME}.edges.tmp.txt"

sed -E "s/.{$bases}/&\n/g" "${NAME}.edges.tmp.txt" | sed '/^$/d' > "${NAME}.edges.txt"
rm "${NAME}.edges.tmp.txt"

echo "Task 3: runBioawk - unknown Ns, headers, lengths"
bioawk -c fastx '{n=gsub(/N/, "", $seq); print $name "\t" n}' "$ASM" > "${NAME}.unknown.txt"
bioawk -c fastx '{print ">"$name"_start";print ">"$name"_end";}' "$ASM" > "${NAME}.headers.txt"
bioawk -c fastx '{print $name, length($seq)}' "$ASM" > "${NAME}.lengths.txt"

echo "Task 4: combineIntoFasta (paste) - headers + edges to FASTA"
paste -d '\n' "${NAME}.headers.txt" "${NAME}.edges.txt" > "${NAME}.combined.fasta"


echo "Task 5: runNCRF - detect TTAGGG telomeres in ends"
echo "Swapping mamba envs..."
mamba deactivate
mamba activate ncrf-py2

NCRF --stats=events TTAGGG < "${NAME}.combined.fasta" | ncrf_summary > "${NAME}.telomeric.summary.txt"
tail -n +2 "${NAME}.telomeric.summary.txt" | cut -f3 | sort | uniq > "${NAME}.telomeric.ends.txt"

mamba deactivate
echo "Returnin to T2T env..."
mamba activate t2t-assess

echo "Task 6: runMashMap - align to CHM13 for chrom assignment (PAF format)"
mashmap --threads 4 --perc_identity 95 --noSplit \
  -r "$REF" -q "$ASM" -o "${NAME}.mashmap.txt"

echo "Task 7: assessCompleteness"

while read -r line; do
  contig=$(echo "$line" | sed 's/_.*//')
  grep "$contig" "${NAME}.lengths.txt" | tr -d '\n'; echo -ne '\t'
  grep "$contig" "${NAME}.unknown.txt" | tr -d '\n'; echo -ne '\t'
  grep "$contig" "${NAME}.mashmap.txt" | tr -s '[:blank:]' '\t'
done < "${NAME}.telomeric.ends.txt" \
  | cut -f1,2,4,10 \
  | sort | uniq -c | sort -rgk1 \
  | tr -s '[:blank:]' '\t' | sed 's/^\t//' \
  > "${NAME}.SUMMARY.txt"

echo -e 'name\tlength\tNs\tchromosome' > "${NAME}.T2T.scaffolds.txt"
echo -e 'name\tlength\tNs\tchromosome' > "${NAME}.T2T.contigs.txt"

if grep -q '^2' "${NAME}.SUMMARY.txt"; then
  grep '^2' "${NAME}.SUMMARY.txt" | awk '{if ($4>=0) print;}' | cut -f2- | sort -k4,4 -V -s >> "${NAME}.T2T.scaffolds.txt" || true
  grep '^2' "${NAME}.SUMMARY.txt" | awk '{if ($4==0) print;}' | cut -f2- | sort -k4,4 -V -s >> "${NAME}.T2T.contigs.txt" || true
else
  echo 'None of the contigs/scaffolds contained telomeres on both ends.' >> "${NAME}.T2T.scaffolds.txt"
fi
echo 'Done.'

echo "Final counts"
CONTIGS=$(tail -n +2 "${NAME}.T2T.contigs.txt" | wc -l)
SCAFFS=$(tail -n +2 "${NAME}.T2T.scaffolds.txt" | wc -l)
echo "T2T contigs: $CONTIGS"
echo "T2T scaffolds: $SCAFFS"
echo "Check ${NAME}.SUMMARY.txt for details (grep '^2' for both-end hits)."
