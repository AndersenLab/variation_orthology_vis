#!/bin/bash

#SBATCH -J pavtrim_indels_invs      
#SBATCH -A eande106         
#SBATCH -p parallel         
#SBATCH -t 2:00:00         
#SBATCH -N 1               
#SBATCH -c 6       
#SBATCH --mail-user=jbauman7@jh.edu  
#SBATCH --mail-type=END     
#SBATCH --output=/vast/eande106/projects/Josh/SV_calling/processed_data/slurm_output/pavtrim_indels.oe  
#SBATCH --error=/vast/eande106/projects/Josh/SV_calling/processed_data/slurm_output/pavtrim_indels.rr  




#Trim pav merged file
#To output in TSV format easily readable in R


bcftools view \
    -m2 -M2 \
    -i '(SVTYPE="INS" || SVTYPE="DEL" || SVTYPE="INV") && 
        FILTER="PASS"' \
    /vast/eande106/projects/Josh/SV_calling/raw_data/c_elegans/pav_vcfs/all_pav_merged_elegans.vcf.gz | awk '
    BEGIN {OFS="\t"} 
    /^##/ {next}
    /^#/ {printf "%s\t%s\t%s\t%s\t%s\t%s", $1, $2, $4, $5, "SV", "SV_length"
            for (i = 10; i <= NF; i++) {
                printf "\t%s", $i
            }
            printf "\n"; next}
    {
        split($3, id, "-");
        sv_type = id[3]
        sv_len = id[4]

        if ((sv_type == "DEL" || sv_type == "INS" || sv_type == "INV" ) && sv_len > 50) {
            printf "%s\t%s\t%s\t%s\t%s\t%s", $1, $2, $4, $5, sv_type, sv_len
            for (i = 10; i <= NF; i++) {
                printf "\t%s", $i
            }
            printf "\n"
        }
    }' >  pav_indels_invs_merged.tsv

