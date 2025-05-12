#!/bin/bash

#SBATCH -J blastcat      
#SBATCH -A eande106         
#SBATCH -p parallel         
#SBATCH -t 06:00:00         
#SBATCH -N 1               
#SBATCH -c 24         
#SBATCH --output=/vast/eande106/projects/Josh/alignment_BLAST/slurm_output/blastcat.oe  
#SBATCH --error=/vast/eande106/projects/Josh/alignment_BLAST/slurm_output/blastcat.rr  



for file in *_wi_blast_N2lib.pb.txt; do
    strain=$(basename $file | cut -d'.' -f1)
    awk -v strain=$strain '{
    print $0 "\t" strain
    }' $file 
done > all_wi_blast_N2_hits.tsv


for file in *_N2_blast_wilib.pb.txt; do
    strain=$(basename $file | cut -d'.' -f1)
    awk -v strain=$strain '{
    print $0 "\t" strain
    }' $file 
done > all_N2_blast_wi_hits.tsv