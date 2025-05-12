#!/bin/bash

#SBATCH -J variantfilter      
#SBATCH -A eande106         
#SBATCH -p parallel         
#SBATCH -t 2:00:00         
#SBATCH -N 1               
#SBATCH -c 6       
#SBATCH --mail-user=jbauman7@jh.edu  
#SBATCH --mail-type=END     
#SBATCH --output=/vast/eande106/projects/Josh/SV_calling/processed_data/slurm_output/variantfilter.oe  
#SBATCH --error=/vast/eande106/projects/Josh/SV_calling/processed_data/slurm_output/variantfilter.rr  



#Arg1: species
#Arg2: chrom
#Arg3: start
#Arg4: stop
#Arg5: strain (character)

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#NOTE: this script is optimized for C. elegans, and needs its functionality to be updated for briggsae and tropicalis. 
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if [[ $1 == "c_elegans" ]]; then
    GATK_vcf="/vast/eande106/projects/Josh/SV_calling/raw_data/c_elegans/WI.20250331.hard-filter.isotype.vcf.gz"
    output_dir="/vast/eande106/projects/Josh/SV_calling/processed_data/vcf_filtering/c_elegans"
    PAV_vcf="/vast/eande106/projects/Josh/SV_calling/raw_data/c_elegans/pav_vcfs/all_pav_merged_elegans.vcf.gz"
    bcsq_snv_annot_vcf="/vast/eande106/projects/Josh/SV_calling/raw_data/c_elegans/CSQ_c_elegans_WBGeneID_GRANTHAM_BLOSUM_test.tsv"
    annovar_snv_annot_vcf="/vast/eande106/projects/Josh/SV_calling/raw_data/c_elegans/ANNOVAR_c_elegans_WBGeneID_GRANTHAM_BLOSUM_test.tsv"
    vep_snv_annot_vcf="/vast/eande106/projects/Josh/SV_calling/raw_data/c_elegans/VEP_c_elegans_WBGeneID_GRANTHAM_BLOSUM_test.tsv"    
elif [[ $1 == "c_troopicalis" ]]; then
    GATK_vcf="/vast/eande106/projects/Josh/SV_calling/raw_data/c_tropicalis/WI.20250331.hard-filter.isotype.vcf.gz"
    output_dir="/vast/eande106/projects/Josh/SV_calling/processed_data/vcf_filtering/c_tropicalis"
elif [[ $1 == "c_briggsae" ]]; then
    GATK_vcf="/vast/eande106/projects/Josh/SV_calling/raw_data/c_briggsae/WI.20250331.hard-filter.isotype.vcf.gz"
    output_dir="/vast/eande106/projects/Josh/SV_calling/processed_data/vcf_filtering/c_briggsae"
else
    echo "Please provide a supported species: c_elegans, c_tropicalis, c_briggsae"
fi




#Filtered VCF 1
#Biallelic snps 
#Did the snp filtering in bcftools
#Next, I used awk to filter just lines with 0/1 or 1/1 variant genotypes (so that only variants in strain of interest are preseent)


bcftools view \
    -r $2:$3-$4 \
    -m2 -M2 -v snps \
    -s $5 \
    $GATK_vcf | awk -v strain="$5" '
    BEGIN {OFS="\t"} 
    /^##/ {next}
    /^#/ {print "CHROM", $2, $4, $5, strain; next}
    {
        split($10, metrics, ":");
        if (metrics[1] != "0/0" && metrics[1] != "./.") {
            print $1, $2, $4, $5, metrics[1]
        }
    }' >  $output_dir/biallelic_snps_${1}_${2}_${3}_${4}_${5}.tsv

#Get annotated snp files

#bscq
grep "$5" $bcsq_snv_annot_vcf | \
    awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5}' | \
    uniq | \
    awk -v chr="$2" -v start="$3" -v end="$4" '$1 == chr && $2 <= end && $2 >= start' > $output_dir/bcsq_annotated_snvs_${1}_${2}_${3}_${4}_${5}.tsv
#annovar
grep "$5" $annovar_snv_annot_vcf | \
    awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $6}' | \
    uniq | \
    awk -v chr="$2" -v start="$3" -v end="$4" '$1 == chr && $2 <= end && $2 >= start' > $output_dir/annovar_annotated_snvs_${1}_${2}_${3}_${4}_${5}.tsv
#vep
grep "$5" $vep_snv_annot_vcf | \
    awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5}' | \
    uniq | \
    awk -v chr="$2" -v start="$3" -v end="$4" '$1 == chr && $2 <= end && $2 >= start' > $output_dir/vep_annotated_snvs_${1}_${2}_${3}_${4}_${5}.tsv



#Filtered VCF 2 (indels < 50bp)
#new header: CHROM, POS, REF, ALT, SV, SV_length, strain genotype
#First use bcftools to subset to just CB4586, then just variations where strain GT is 0/1 or 1/1
#Within these lines, define the length of the reference and alt alleles 
#Skip indels greater than 50 
#assign deletion or insertion based on length of reference and alt alleles 
#Print the new line 

bcftools view \
    -r $2:$3-$4 \
    -m2 -M2 -v indels \
    -s $5 \
    $GATK_vcf | awk -v strain="$5" '
    BEGIN {OFS="\t"} 
    /^##/ {next} 
    /^#/ {print "CHROM", $2, $4, $5, "Variant", "Var_length", strain; next}
    {
        split($10, metrics, ":");
        if (metrics[1] != "0/0" && metrics[1] != "./.") {

            ref_len=length($4); alt_len=length($5);

            if (ref_len >= 50 || alt_len >= 50) {next;}

            if (ref_len > alt_len) {
                sv_type = "DEL";
                sv_len = ref_len - alt_len;
            } else if (ref_len < alt_len) {
                sv_type = "INS";
                sv_len = alt_len - ref_len;
            } else {
                sv_type = "SNP"; 
                sv_len = 0;
            }

            print $1, $2, $4, $5, sv_type, sv_len, metrics[1]
        }
    }' > $output_dir/biallelic_indels_under50bp_${1}_${2}_${3}_${4}_${5}.tsv


##Now, get PAV SV's using the merged file I created earlier
#Select region of interest, biallelic variants, in strain of interest, and then filter for only SVs that pass quality check
#Finall, create a tsv of just indels > 50bp. 

bcftools view \
    -r $2:$3-$4 \
    -m2 -M2 \
    -s $5 \
    -i '(SVTYPE="INS" || SVTYPE="DEL" || SVTYPE="INV" || SVTYPE="TRA") && 
        FILTER="PASS"' \
    $PAV_vcf | awk -v strain="$5" '
    BEGIN {OFS="\t"} 
    /^##/ {next}
    /^#/ {print "CHROM", $2, $4, $5, "SV", "SV_length", strain; next}
    {
        split($3, id, "-");
        sv_type = id[3]
        sv_len = id[4]

        if ($10 != "./." && (sv_type == "DEL" || sv_type == "INS" || sv_type == "INV") && sv_len > 50) {
            print $1, $2, $4, $5, sv_type, sv_len, $10
        }
    }' >  $output_dir/pav_SVs_${1}_${2}_${3}_${4}_${5}.tsv


#Put all this information into the Rscript for variant plotting 
Rscript /vast/eande106/projects/Josh/SV_calling/scripts/variantvis_josh_strainspecific.R $1 $2 $3 $4 $5
