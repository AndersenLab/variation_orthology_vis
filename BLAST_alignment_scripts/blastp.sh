#!/bin/bash

#SBATCH -J blastp      
#SBATCH -A eande106         
#SBATCH -p parallel         
#SBATCH -t 06:00:00         
#SBATCH -N 1               
#SBATCH -c 24         
#SBATCH --output=/vast/eande106/projects/Josh/alignment_BLAST/slurm_output/blastp.oe  
#SBATCH --error=/vast/eande106/projects/Josh/alignment_BLAST/slurm_output/blastp.rr  

source activate blast

#WI Protein against N2 library
blastp -query $prot -db /vast/eande106/projects/Josh/alignment_BLAST/reference/libraries/protein_lib/N2.WBonly.WS283.PConly.prot.fa -out /vast/eande106/projects/Josh/alignment_BLAST/blastp/${prot%.*}_wi_blast_N2lib.pb.txt -outfmt 6 -num_threads 24

#N2 protein against WI library 
blastp -query /vast/eande106/projects/Josh/alignment_BLAST/reference/protein/N2.WBonly.WS283.PConly.prot.fa -db /vast/eande106/projects/Josh/alignment_BLAST/wi/libraries/${prot} -out /vast/eande106/projects/Josh/alignment_BLAST/blastp/${prot%.*}_N2_blast_wilib.pb.txt -outfmt 6 -num_threads 24

#cd into /vast/eande106/projects/Josh/alignment_BLAST/wi/protein directory and run: 
# for file in *protein.fa; do sbatch --export=prot=$file /vast/eande106/projects/Josh/alignment_BLAST/scripts/blastp.sh; done


