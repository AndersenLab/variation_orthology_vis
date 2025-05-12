#!/bin/bash

#SBATCH -J blastx      
#SBATCH -A eande106         
#SBATCH -p parallel         
#SBATCH -t 06:00:00         
#SBATCH -N 1               
#SBATCH -c 24         
#SBATCH --output=/vast/eande106/projects/Josh/alignment_BLAST/slurm_output/blastx.oe  
#SBATCH --error=/vast/eande106/projects/Josh/alignment_BLAST/slurm_output/blastx.rr  

source activate blast

#WI Protein against N2 library
blastx -query $prot -db /vast/eande106/projects/Josh/alignment_BLAST/reference/libraries/genome_lib/c_elegans.PRJNA13758.WS283.genome.fa -out /vast/eande106/projects/Josh/alignment_BLAST/blastx/${prot%.*}_wi_blast_N2genlib.blastx.txt -outfmt 6 -num_threads 24

#cd into /vast/eande106/projects/Josh/alignment_BLAST/wi/protein directory and run: 
# for file in *protein.fa; do sbatch --export=prot=$file /vast/eande106/projects/Josh/alignment_BLAST/scripts/blastx.sh; done