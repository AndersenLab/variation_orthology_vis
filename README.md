# variation_orthology_vis
Josh Bauman rotation spring 2025 - wrote bash scripts to filter GATK and PAV VCFs, R scripts for visualizing this variation, and worked on validating BRAKER gene model calls between N2 and WS using BLASTP 

## Structural variant calling

### `pav_hdr_intersect.sh` 
This script is necessary as part of the `PAV_INDEL_vis_josh.R` script to prune out HDRs from the SV plot.

### `PAV_INDEL_vis_josh.R`  
This script plots all the PAV results from long read assemblies (from `pavmergetrim.sh`) faceted by chromosome. HDRs are pruned with `pav_hdr_intersect.sh`

### `pavmerge2.sh`  
Merges all PAV results when run in the same directory that contains all of the PAV output vcfs

### `pavmerge_trim_invs.sh`  
Same as `pavmerge_trim.sh`, but includes inversions

### `pavmerge_trim.sh`  
converts the output of pavmerge2 into a format that can be read by `PAV_INDEL_vis_josh.R`

### `variantvis_josh_strainspecific.R` 
Rscript that takes the output of `vcf_filtering` and creates the variant plot, annotated SNV plot, and data tables. It is run internally by `vcf_filtering.sh`

### `vcf_filtering.sh` 
the main function to run to generate plot and variant data tables. Arguments are made in this format (glc-1 locus in this case): `sbatch vcf_filtering.sh c_elegans V 16115967 16276907 CB4856`

### `variantreport.Rmd`
This script takes the `.Rdata` file from `variantvis_josh_strainspecific` and creates a markdown report with the variant plot and variant tables. This was done as a workaround because I was unable to knit an html on the Rstudio server. Ideally this script would be revised and all the code from `variantvis_josh_strainspecific.R` would be run in the first markdown cell. That would require this markdown doc to inherit arguments from the `vcf_filtering.sh` script, which might be tricky. 


## Blast alignment annotations

### `blastcat.sh`
Concatenates the BLASTp results for all the wi->N2 protein blasts and the N2->wi protein blasts, and adds a column for strain

### `blastp.sh`
Performs the blast search for each wi protein fasta against the N2 database (and the N2 fasta agains each wi database). Results go into `blastcat.sh`

### `blastx.sh`  
performs BLASTx against each wi protein fasta against the N2 genome 

### `prot_libs.sh`  
Makes blast databases for each wi protein fasta 

### `RBBH.R`
A rudimentary script for performing RBBH on CB4856 protein fast -> N2 protein fasta and reciprocal blast. Needs to be modularized to work on any strain, and needs to be updated to account for cases where the RBBH proteins are not always the highest hits sorted by ascending Evalue and descending bitscore. 