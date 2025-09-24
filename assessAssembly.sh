### Re-implemention of https://github.com/biomonika/HPP/blob/main/assembly/wdl/workflows/assessAsemblyCompletness.wdl outside of Cromwell

!/bin/bash


set -euo pipefail  # Fail on errors

ASM="hap2-dualscaff.fa"      # Swap to hap2.fa for second run
REF="chm13v2.0.fa"
NAME=$(basename "$ASM" .fa)  # e.g., hap1

echo "Running T2T assessment for $ASM..."

# Task 1: formatAssembly (seqtk) - rewrap to single-line FASTA
docker run --user 1004:1004 --rm -v $(pwd):/data -w /data quay.io/biocontainers/seqtk:1.3--hed695b0_2 \
  seqtk seq -l0 /data/"$ASM" > "${NAME}".formatted.fa || { echo "Task 1 failed"; exit 1; }

# Task 2: generateAssemblyEdges (gawk) - extract 1000 bp ends
docker run --user 1004:1004 --rm -v $(pwd):/data -w /data quay.io/biocontainers/gawk:5.1.0--2 bash -c "
  set -euo pipefail
  bases=1000
  cat /data/'${NAME}'.formatted.fa | gawk -v len=\"\$bases\" -F '' '{if (NR % 2 == 0) {for (i=1; i<=NF; i++) {if (i <= len || (NF-i) < len) {printf \$(i)}}; printf \"\\n\"}}' > /data/'${NAME}'.edges.tmp.txt
  cat /data/'${NAME}'.edges.tmp.txt | sed -e \"s/.{\${bases}}/&\\n/g\" | sed '/^$/d' > /data/'${NAME}'.edges.txt
  rm /data/'${NAME}'.edges.tmp.txt
" || { echo "Task 2 failed"; exit 1; }

# Task 3: runBioawk - unknown Ns, headers, lengths
docker run --user 1004:1004 --rm -v $(pwd):/data -w /data quay.io/biocontainers/bioawk:1.0--h577a1d6_13 bash -c "
  set -euo pipefail
  bioawk -c fastx '{n=gsub(/N/, \"\", \$seq); print \$name \"\\t\" n}' /data/'$ASM' > /data/'${NAME}'.unknown.txt
  bioawk -c fastx '{print \">\"\$name\"_start\";print \">\"\$name\"_end\";}' /data/'$ASM' > /data/'${NAME}'.headers.txt
  bioawk -c fastx '{ print \$name, length(\$seq) }' < /data/'$ASM' > /data/'${NAME}'.lengths.txt
" || { echo "Task 3 failed"; exit 1; }

# Task 4: combineIntoFasta (paste) - headers + edges to FASTA
docker run --user 1004:1004 --rm -v $(pwd):/data -w /data ubuntu \
  bash -c "
    set -euo pipefail
    paste -d '\n' /data/'${NAME}'.headers.txt /data/'${NAME}'.edges.txt > /data/'${NAME}'.combined.fasta
  " || { echo "Task 4 failed"; exit 1; }

# Task 5: runNCRF - detect TTAGGG telomeres in ends
docker run --user 1004:1004 --rm -v $(pwd):/data -w /data --memory=16g quay.io/biocontainers/ncrf:1.01.02--h7b50bb2_6 bash -c "
  set -euo pipefail
  cat /data/'${NAME}'.combined.fasta | NCRF --stats=events TTAGGG | ncrf_summary > /data/'${NAME}'.telomeric.summary.txt
  tail -n +2 /data/'${NAME}'.telomeric.summary.txt | cut -f3 | sort | uniq > /data/'${NAME}'.telomeric.ends.txt
" || { echo "Task 5 failed"; exit 1; }

# Task 6: runMashMap - align to CHM13 for chrom assignment (PAF format)
docker run --user 1004:1004 --rm -v $(pwd):/data -w /data --memory=24g --cpus=4 quay.io/biocontainers/mashmap:3.1.3--pl5321hb4818e0_2 \
  mashmap --threads 4 --perc_identity 95 --noSplit -r /data/"$REF" -q /data/"$ASM" -o /data/"${NAME}".mashmap.txt || { echo "Task 6 failed"; exit 1; }

# Task 7: assessCompletness - build SUMMARY, filter T2T (both ends=2, contigs: Ns=0, scaffolds: Ns>0; chrom from mashmap f10)
docker run --user 1004:1004 --rm -v $(pwd):/data -w /data ubuntu bash -c "
  set -euo pipefail
  cat /data/'${NAME}'.telomeric.ends.txt | while read line; do
    contig=\$(echo \$line | sed 's/_.*//')
    egrep \"\$contig\" /data/'${NAME}'.lengths.txt | tr -d '\\n'; echo -ne '\\t'
    egrep \"\$contig\" /data/'${NAME}'.unknown.txt | tr -d '\\n'; echo -ne '\\t'
    egrep \"\$contig\" /data/'${NAME}'.mashmap.txt | tr -s '[:blank:]' '\\t'
  done | cut -f1,2,4,10 | sort | uniq -c | sort -rgk1 | tr -s '[:blank:]' '\\t' | sed 's/^\\t//' > /data/'${NAME}'.SUMMARY.txt

  echo -e 'name\\tlength\\tNs\\tchromosome' > /data/'${NAME}'.T2T.scaffolds.txt
  echo -e 'name\\tlength\\tNs\\tchromosome' > /data/'${NAME}'.T2T.contigs.txt

  if grep -q '^2' /data/'${NAME}'.SUMMARY.txt; then
    cat /data/'${NAME}'.SUMMARY.txt | grep '^2' | awk '{if (\$4>=0) print;}' | cut -f2- | sort -k4,4 -V -s >> /data/'${NAME}'.T2T.scaffolds.txt || true
    cat /data/'${NAME}'.SUMMARY.txt | grep '^2' | awk '{if (\$4==0) print;}' | cut -f2- | sort -k4,4 -V -s >> /data/'${NAME}'.T2T.contigs.txt || true
  else
    echo 'None of the contigs/scaffolds contained telomeres on both ends.' >> /data/'${NAME}'.T2T.scaffolds.txt
  fi
  echo 'Done.'
" || { echo "Task 7 failed"; exit 1; }

# Final counts (match preprint Table 1)
CONTIGS=$(tail -n +2 "${NAME}".T2T.contigs.txt | wc -l)
SCAFFS=$(tail -n +2 "${NAME}".T2T.scaffolds.txt | wc -l)
echo "T2T contigs: $CONTIGS"
echo "T2T scaffolds: $SCAFFS"
echo "Check ${NAME}.SUMMARY.txt for details (grep '^2' for both-end hits)."
