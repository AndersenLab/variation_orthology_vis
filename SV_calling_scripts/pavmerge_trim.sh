#!/bin/bash
#Trim pav merged file
#For plotting all indels, faceted by chromosome, in R


bcftools view \
    -m2 -M2 \
    -i '(SVTYPE="INS" || SVTYPE="DEL") && 
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

        if ((sv_type == "DEL" || sv_type == "INS") && sv_len > 50) {
            printf "%s\t%s\t%s\t%s\t%s\t%s", $1, $2, $4, $5, sv_type, sv_len
            for (i = 10; i <= NF; i++) {
                printf "\t%s", $i
            }
            printf "\n"
        }
    }' >  pav_indels_merged.tsv

