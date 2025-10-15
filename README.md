## T2T Assm + QC
Usage:  

``./createEnvs.sh`` # Set up required Mamba environments

Requires Hs1 reference in the script root directory as "hs1.fa".  
`` wget https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz ``  
`` gunzip hs1.fa.gz ``

``./t2t_assess <haplotype_fa>`` # haplotype fasta file, can be produced via awk (Hifiasm FAQs).
