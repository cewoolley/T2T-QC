## T2T Assm + QC
Usage:  

``./createEnvs.sh`` # Set up required Mamba environments

Requires Hs1 reference in the script root directory as "hs1.fa".  
`` wget https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz ``  
`` gunzip hs1.fa.gz ``

``./t2t_assess <haplotype_fa>`` # haplotype fasta file, can be produced via awk (Hifiasm FAQs).  

``./t2t_batch_assess.sh -r <analysis_dir, optional> -R <hs1.fa path, optional> -t <threads, defaults to 8> `` # Scaled up version of the above, iterates across haplotype contig fasta files and QCsinto a summary table.
