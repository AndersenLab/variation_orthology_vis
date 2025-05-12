#!/bin/bash


#SBATCH -J pavmerge    
#SBATCH -A eande106         
#SBATCH -p parallel         
#SBATCH -t 2:00:00         
#SBATCH -N 1               
#SBATCH -c 6       
#SBATCH --mail-user=jbauman7@jh.edu  
#SBATCH --mail-type=END     
#SBATCH --output=/vast/eande106/projects/Josh/SV_calling/processed_data/slurm_output/pavmerge.oe  
#SBATCH --error=/vast/eande106/projects/Josh/SV_calling/processed_data/slurm_output/pavmerge.rr  




#CD into directory with all PAV output vcfs and run this script
 
#convert genotype field from "1" to "1/1"

for file in *.vcf.gz; do
    name=${file%.vcf.gz}
    bcftools view $file | awk '
    BEGIN {OFS="\t"} 
    /^##/ {print $0; next}
    /^#/ {print $0; next}
    {
        if($10 == 1){
        genotype = "1/1"
        } else {
        genotype = "./."
        }
        print $1, $2, $3, $4, $5, $6, $7, $8, $9, genotype
    }' > gt_corrected_${name}.vcf

    #Compress the file
    bcftools view -Oz -o gt_corrected_${name}.vcf.gz gt_corrected_${name}.vcf
    
    #generate index for each file
    bcftools index gt_corrected_${name}.vcf.gz

    echo "$file done"
done 

#Next, merge all filtered vcfs, generate index, and clean directory
bcftools merge -o all_pav_merged_elegans.vcf.gz gt_corrected_*.vcf.gz 
bcftools index all_pav_merged_elegans.vcf.gz

rm gt_corrected_*.vcf*

