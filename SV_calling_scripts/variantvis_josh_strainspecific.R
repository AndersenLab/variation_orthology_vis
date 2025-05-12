library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(cowplot)
library(ape)
library(grid)
library(stringr)
library(ggthemes)

##Make sure to clear objects between runs of this script in Rstudio. If running through terminal, should be fine. 

#Run commands below if you just want to test this script without running the vcf_filtering bash script
#args <- c("c_elegans", "V", "16115967", "16276907", "CB4856")
#args <- c("c_elegans", "V", "16185202", "16185502", "CB4856")


##JOSH##
#Inherit arguments from bash script 
args <- commandArgs(trailingOnly = TRUE)
args_species <- args[1] 
args_chrom <- args[2]
args_start <- args[3]
args_stop <- args[4]
args_strain <- args[5]

#Find the length of the interval (to keep polygon widths consistent when plotting)
interval_length <- as.numeric(args_stop) - as.numeric(args_start) 
#Set a size buffer for plotting very small variants
size_buffer <- 0.0001 * interval_length


##### This section of code is taken from Nic's old haplotypePlotter R script ####

#read collapsed reference HDRs
collapsed_ff <- readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/glc1_variation/HDR_haplotypePlotter/input/HDR_5kbclust_collapsed_wFreq.tsv")

#read strain-specific HDRs
# all_SR_calls <- readr::read_tsv("./input/HDR_allStrain_5kbclust_1IBfilt.tsv")

#read all pairwise genome coordinate comparisons
### This will need to eventually change to all pairwise alignments among all WSs to the reference
transformed_coords <- readr::read_tsv("/vast/eande106/projects/Nicolas/c.elegans/reference_genealn/N2vCB/N2_hifi_transformed2.tsv",col_names = F) %>% dplyr::mutate(STRAIN="CB4856")

colnames(transformed_coords) <- c("S1","E1","S2","E2","L1","L2","IDY","LENR","LENQ","REF","HIFI","STRAIN")

#read concatentated gene models of every genome
gffCat <- ape::read.gff("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/glc1_variation/HDR_haplotypePlotter/input/all_LRiso.gff")
gffCat1 <- ape::read.gff("/vast/eande106/projects/Nicolas/WI_PacBio_genomes/annotation/elegans/braker_runs/gff/CB4856.braker.gff3") %>% dplyr::mutate(STRAIN="CB4856")
#gffCat2 <- ape::read.gff("/vast/eande106/projects/Nicolas/c.elegans/N2/wormbase/WS283/N2.WBonly.WS283.PConly.gff3") %>% dplyr::mutate(STRAIN="N2")
gffCat2 <- ape::read.gff("/vast/eande106/data/c_elegans/genomes/PRJNA13758/WS283/csq/c_elegans.PRJNA13758.WS283.csq.gff3") %>% dplyr::mutate(STRAIN="N2")

gffCat <- rbind(gffCat1,gffCat2)


#set your target  coordinates
#GLC-1
hdr_chrom = args_chrom
hdr_start_pos = as.numeric(args_start)
hdr_end_pos = as.numeric(args_stop)

#offset lets you explore adjacent regions
offset = 0
hap_chrom = hdr_chrom
hap_start = hdr_start_pos #- offset
hap_end = hdr_end_pos #+ offset 

#use reference coordinates from g2g alginments to pull the contigs that contain the alt haplotypes for the HDR
# hap_coords <- transformed_coords %>%
#   dplyr::filter((REF == hap_chrom & hap_start >= S1 & hap_start <= E1 ) | 
#                   (REF == hap_chrom & hap_end >= S1 & hap_end <= E1) | 
#                   (REF == hap_chrom & S1 >= hap_start & E1 <= hap_end)) %>%
#   dplyr::mutate(inv=ifelse(S2>E2,T,F)) 
# # dplyr::select(-S2,-E2) %>%
# # dplyr::rename(S2=newS2,E2=newE2)
# 
# 
# #keep only the contig with the largest extent of alignment with the REF HDR
# tigFilt <- hap_coords %>% 
#   dplyr::group_by(STRAIN) %>%
#   dplyr::mutate(nalign = n()) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(HIFI) %>%
#   dplyr::mutate(ntig= n()) %>%
#   dplyr::mutate(tigsize=sum(L2)) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(STRAIN) %>%
#   dplyr::filter(tigsize == max(tigsize)) %>%
#   dplyr::ungroup()
# 
# 
# #get the minimum and maximum boundary of the WILD genome alignments that contain the HDR
# HV_boundary <- tigFilt %>%
#   dplyr::mutate(refStart=min(S1),refEnd=max(E1)) %>%
#   dplyr::group_by(STRAIN) %>%
#   dplyr::mutate(boundStart=min(S2), boundEnd=max(E2)) %>%
#   dplyr::distinct(STRAIN, .keep_all = T) %>%
#   dplyr::select(HIFI,boundStart,boundEnd,STRAIN,REF,refStart,refEnd) %>%
#   dplyr::ungroup()
# 
# #filter the concatenated GFF to extract the gene models of each WILD genome contig boundary
# filtGff <- gffCat %>%
#   dplyr::filter(type=="gene") %>%
#   dplyr::filter(seqid %in% HV_boundary$HIFI) %>%
#   dplyr::left_join(HV_boundary,by=c('seqid'='HIFI',"STRAIN")) %>%
#   dplyr::filter((start >= boundStart & start <= boundEnd) | 
#                   (end >= boundStart & end <= boundEnd)) %>%
#   dplyr::group_by(STRAIN) %>%
#   dplyr::mutate(ngene=n()) %>%
#   dplyr::arrange(ngene) %>%
#   dplyr::mutate(gid=cur_group_id()) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(boundSize=abs(boundStart-boundEnd))

#extract the REF genes and pseudogenes
N2Start = hap_start
N2End = hap_end
N2_genes <- gffCat %>%
  dplyr::filter(seqid==hap_chrom & type=="gene") %>%
  dplyr::mutate(refStart=N2Start,refEnd=N2End) %>%
  dplyr::filter((start >= hap_start & start <= hap_end) | 
                  (end >= hap_start & end <= hap_end)) %>%
  dplyr::filter(grepl("biotype=protein_coding|biotype=polymorphic_pseudogene",attributes)) %>%
  tidyr::separate(attributes,into=c("pre","post"),sep=';sequence_name=') %>%
  tidyr::separate(post,into=c("seqname","post2"),sep=';biotype=') %>%
  dplyr::mutate(seqname=paste0("Transcript_",seqname)) %>%
  tidyr::separate(pre,into=c("ID","Name","rest2"),sep=";") %>%
  dplyr::mutate(Name=gsub("Name=","",Name)) %>%
  dplyr::mutate(
    genetype = ifelse(grepl("polymorphic_pseudogene", post2), "pseudogene", "protein_coding")
  ) 



#get N2 genes
#Lance seemed to be doing some shifting of the start and stop. I changed to just use raw start and stop, for keep axes consistent in plotting. 
N2ad <- N2_genes %>%
  dplyr::mutate(STRAIN="N2",ortho_status=T) %>%
  dplyr::rename(N2=seqname) %>%
  #dplyr::mutate(shift=min(start)) %>% #nic
  #dplyr::mutate(newstart=start-shift,newend=end-shift) %>% #nic
  dplyr::mutate(newstart=start,newend=end) %>%
  dplyr::select(seqid,STRAIN,N2,ortho_status,newstart,newend,genetype,strand) %>%
  dplyr::mutate(sp = ifelse(strand == "+", 2.25,1.75)) 









######## Here is where the TSVs you create in bash will be directly loaded in and used ######## 

# Load in SV and SNV calls by paftools (Not needed right now)


# df = readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/elegans/vcf/elegans.merged.1kbCOV.5kbALIGN.annotated.final.vcf")

# CB4856var <- df %>%
#   dplyr::rename(CHROM = '#CHROM', start = POS) %>%
#   dplyr::select(CHROM,start,end,REF,ALT,INFO,CB4856) %>%
#   dplyr::filter(CHROM == "V", (start >= 16115967 & start <= 16276907)) %>%
#   dplyr::filter(CB4856 != "./.")
# 
# plot_df <- CB4856var %>%
#   dplyr::mutate(
#     lenDEL = case_when(INFO == "DEL" ~ (end - start), TRUE ~ NA_real_),
#     lenINS = case_when(INFO == "INS" ~ (end - start), TRUE ~ NA_real_))

# SNPs <- plot_df %>%
#   dplyr::filter(INFO == "SNP") %>%
#   dplyr::select(CHROM,start,end,REF,ALT,INFO,CB4856) %>%
#   dplyr::filter(!grepl(",", ALT))

# deletions_paf <- plot_df %>%
#   dplyr::filter(INFO == "DEL") %>%
#   dplyr::select(CHROM,start,end,REF,ALT,INFO,CB4856,lenDEL) %>%
#   dplyr::filter(lenDEL >= 50)
# 
# insertions_paf <- plot_df %>%
#   dplyr::filter(INFO == "INS") %>%
#   dplyr::select(CHROM,start,end,REF,ALT,INFO,CB4856,lenINS) %>%
#   dplyr::filter(lenINS >= 50)



# Load in SNPs called by GATK pipeline with SR data
#SR = readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/elegans/vcf/SR.GATK.CB4856.final.vcf")

##LANCE##
# GATK <- SR %>%
#   dplyr::filter(POS >= 16115967 & POS <= 16276907) %>%
#   dplyr::rename(start = POS) %>%
#   dplyr::mutate(end = start+1) %>%
#   dplyr::filter(CB4856 != './.' & CB4856 != '0/0')

##JOSH##
#Loading in SNP calls from strain of interest
snp_vcf <- paste("/vast/eande106/projects/Josh/SV_calling/processed_data/vcf_filtering", 
                 args_species, paste(
                   "biallelic_snps",
                   "_",
                   args_species, 
                   "_",
                   args_chrom, 
                   "_",
                   args_start, 
                   "_",
                   args_stop, 
                   "_",
                   args_strain,
                   ".tsv",
                   sep = ""
                 ), sep = "/"
)

SR <- readr::read_tsv(snp_vcf)
GATK <- SR  %>%
  dplyr::rename(start = POS) %>%
  dplyr::mutate(end = as.numeric(start)+1)

##Add in SNV variant annotations 
bcsq_annot_snp_vcf <- paste("/vast/eande106/projects/Josh/SV_calling/processed_data/vcf_filtering", 
                 args_species, paste(
                   "bcsq_annotated_snvs",
                   "_",
                   args_species, 
                   "_",
                   args_chrom, 
                   "_",
                   args_start, 
                   "_",
                   args_stop, 
                   "_",
                   args_strain,
                   ".tsv",
                   sep = ""
                 ), sep = "/"
)

bcsq_annot_df <- readr::read_tsv(bcsq_annot_snp_vcf, col_names = FALSE)
colnames(bcsq_annot_df) <- c("chrom", "pos", "ref", "alt", "bcsq_annot")

annovar_annot_snp_vcf <- paste("/vast/eande106/projects/Josh/SV_calling/processed_data/vcf_filtering", 
                            args_species, paste(
                              "annovar_annotated_snvs",
                              "_",
                              args_species, 
                              "_",
                              args_chrom, 
                              "_",
                              args_start, 
                              "_",
                              args_stop, 
                              "_",
                              args_strain,
                              ".tsv",
                              sep = ""
                            ), sep = "/"
)

annovar_annot_df <- readr::read_tsv(annovar_annot_snp_vcf, col_names = FALSE)
colnames(annovar_annot_df) <- c("chrom", "pos", "ref", "alt", "annovar_annot")

vep_annot_snp_vcf <- paste("/vast/eande106/projects/Josh/SV_calling/processed_data/vcf_filtering", 
                            args_species, paste(
                              "vep_annotated_snvs",
                              "_",
                              args_species, 
                              "_",
                              args_chrom, 
                              "_",
                              args_start, 
                              "_",
                              args_stop, 
                              "_",
                              args_strain,
                              ".tsv",
                              sep = ""
                            ), sep = "/"
)

vep_annot_df <- readr::read_tsv(vep_annot_snp_vcf, col_names = FALSE)
colnames(vep_annot_df) <- c("chrom", "pos", "ref", "alt", "vep_annot")

#fix vep dataframe to remove rows with multiple annotations (only keep the first one) 
vep_annot_df <- vep_annot_df %>%
  group_by(chrom, pos, ref, alt) %>%
  slice(1) %>%  # keep just the first annotation
  ungroup()


#Merge the snv dataframe with annotation dataframes

if(nrow(bcsq_annot_df > 0)) {
  GATK_annot <- GATK %>%
    #bcsq annotations
    dplyr::left_join(bcsq_annot_df, by = c("start" = "pos")) %>%
    dplyr::mutate(bcsq_annot = if_else(stringr::str_starts(bcsq_annot, "@"), NA_character_, bcsq_annot)) %>%
    dplyr::mutate(bcsq_annot = if_else(bcsq_annot == "N/A", NA_character_, bcsq_annot)) %>%
    dplyr::mutate(bcsq_annot = if_else(bcsq_annot == "*synonymous", "synonymous", bcsq_annot)) %>%
    dplyr::mutate(bcsq_annot = if_else(bcsq_annot == "*missense", "missense", bcsq_annot)) %>%
    dplyr::mutate(bcsq_annot = if_else(bcsq_annot == "splice_region", "splice region", bcsq_annot)) %>%
    dplyr::mutate(bcsq_annot = if_else(bcsq_annot == "5_prime_utr", "5 prime UTR", bcsq_annot)) %>%
    dplyr::mutate(bcsq_annot = if_else(bcsq_annot == "3_prime_utr", "3 prime UTR", bcsq_annot)) %>%
    dplyr::mutate(bcsq_annot = if_else(bcsq_annot == "non_coding", "noncoding", bcsq_annot)) %>%
    dplyr::mutate(bcsq_annot = if_else(bcsq_annot == "stop_gained", "stop gained", bcsq_annot)) %>%
    dplyr::select(c(CHROM, start, end, REF, ALT, bcsq_annot)) %>%
    #annovar annotations
    dplyr::left_join(annovar_annot_df, by = c("start" = "pos")) %>%
    dplyr::select(c(CHROM, start, end, REF, ALT, bcsq_annot, annovar_annot)) %>%
    dplyr::mutate(annovar_annot = if_else(annovar_annot == "synonymous_SNV", "synonymous", annovar_annot)) %>%
    dplyr::mutate(annovar_annot = if_else(annovar_annot == "nonsynonymous_SNV", "nonsynonymous", annovar_annot)) %>%
    dplyr::mutate(annovar_annot = if_else(annovar_annot == "stopgain", "stop gained", annovar_annot)) %>%
    dplyr::mutate(annovar_annot = if_else(annovar_annot == "N/A", NA_character_, annovar_annot)) %>%
    #VEP annotations
    dplyr::left_join(vep_annot_df, by = c("start" = "pos")) %>%
    dplyr::select(c(CHROM, start, end, REF, ALT, bcsq_annot, annovar_annot, vep_annot))
}
  




#Now, load in indel cells from strain of interest 

##LANCE##
# indels_paf = readr::read_tsv("/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/paftools/elegans/vcf/SR.GATK.CB4856.indels.final.vcf")
# 
# G_del_paf <- indels_paf %>%
#   dplyr::filter(START >= 16115967 & START <= 16276907) %>%
#   dplyr::filter(ANNO == "DEL") %>%
#   dplyr::filter() %>%
#   dplyr::filter((END - START) < 50)
# 
# G_ins_paf <- indels_paf %>%
#   dplyr::filter(START >= 16115967 & START <= 16276907) %>%
#   dplyr::filter(ANNO == "INS") %>%
#   dplyr::filter((END - START) < 50)

##JOSH##
indel_vcf <- paste("/vast/eande106/projects/Josh/SV_calling/processed_data/vcf_filtering", 
                 args_species, paste(
                   "biallelic_indels_under50bp",
                   "_",
                   args_species, 
                   "_",
                   args_chrom, 
                   "_",
                   args_start, 
                   "_",
                   args_stop, 
                   "_",
                   args_strain,
                   ".tsv",
                   sep = ""
                 ), sep = "/"
)
indels <- readr::read_tsv(indel_vcf)

G_del <- indels %>%
  dplyr::filter(Variant == "DEL") %>% 
  dplyr::mutate(START = POS) %>%
  dplyr::mutate(END = as.numeric(POS) + as.numeric(Var_length))

G_ins <- indels %>%
  dplyr::filter(Variant == "INS") %>% 
  dplyr::mutate(START = POS)

#Update insertions with info for plotting correctly 
make_ins_polygon <- function(start, y_pos, triangle_width){
  tibble(
    x = c(
      as.numeric(start)-size_buffer,
      as.numeric(start)+size_buffer,
      as.numeric(start)+size_buffer,
      as.numeric(start) + triangle_width/2,
      as.numeric(start) - triangle_width/2,
      as.numeric(start)-size_buffer
    ),
    y = c(
      y_pos,
      y_pos, 
      y_pos + 0.25,
      y_pos + 0.25 + 0.02, 
      y_pos + 0.25 + 0.02,
      y_pos + 0.25
    )
  )
}

#Run the next line for cases where no insertions are detected (prevents errors with plotting)
if(nrow(G_ins) == 0){
  G_ins_check <- G_ins
  G_ins <- G_ins %>% mutate_all(as.numeric)
  G_ins[1, ] <- as.list(rep(0, ncol(G_ins)))
}

gatk_ins_plotting <- G_ins %>%
  dplyr::rowwise() %>%
  dplyr::mutate(polygon = list(make_ins_polygon(start = as.numeric(START), y_pos = 0.25, triangle_width = 0.005 * interval_length))) %>%
  unnest(polygon)

#Run the next line for cases where no insertions are detected (prevents errors with plotting)
if(exists("G_ins_check")){
  if(nrow(G_ins_check) == 0){
    gatk_ins_plotting[,9] <- 0
    gatk_ins_plotting[,10] <- 0
  }
}

##JOSH##
##Load in INDEL/INV > 50bp calls from PAV
sv_vcf <- paste("/vast/eande106/projects/Josh/SV_calling/processed_data/vcf_filtering", 
                   args_species, paste(
                     "pav_SVs",
                     "_",
                     args_species, 
                     "_",
                     args_chrom, 
                     "_",
                     args_start, 
                     "_",
                     args_stop, 
                     "_",
                     args_strain,
                     ".tsv",
                     sep = ""
                   ), sep = "/"
)

svs <- readr::read_tsv(sv_vcf)

deletions <- svs %>%
  dplyr::filter(SV == "DEL") %>% 
  dplyr::mutate(START = as.numeric(POS)) %>%
  dplyr::mutate(END = as.numeric(POS) + as.numeric(SV_length))

insertions <- svs %>%
  dplyr::filter(SV == "INS") %>% 
  dplyr::mutate(START = as.numeric(POS))

inversions <- svs %>%
  dplyr::filter(SV == "INV") %>% 
  dplyr::mutate(START = as.numeric(POS)) %>%
  dplyr::mutate(END = as.numeric(POS) + as.numeric(SV_length))

#Update insertion dataframe with info for plotting correctly 

make_big_ins_polygon <- function(start, y_pos, triangle_width){
  tibble(
    x = c(
      as.numeric(start)-size_buffer,
      as.numeric(start)+size_buffer,
      as.numeric(start)+size_buffer,
      as.numeric(start) + triangle_width/2,
      as.numeric(start) - triangle_width/2,
      as.numeric(start)-size_buffer
    ),
    y = c(
      y_pos,
      y_pos, 
      y_pos + 0.25,
      y_pos + 0.25 + 0.06, 
      y_pos + 0.25 + 0.06,
      y_pos + 0.25
    )
  )
}

#Run the next line for cases where no insertions are detected (Prevents errors with plotting)
if(nrow(insertions) == 0){
 insertions_check <- insertions
 insertions <- insertions %>% mutate_all(as.numeric)
 insertions[1, ] <- as.list(rep(0, ncol(insertions)))
}

pav_ins_plotting <- insertions %>%
  dplyr::rowwise() %>%
  dplyr::mutate(polygon = list(make_big_ins_polygon(start = as.numeric(START), y_pos = 0.25, triangle_width = 0.01 * interval_length))) %>%
  unnest(polygon)

#Run the next line for cases where no insertions are detected (prevents errors with plotting)
if(exists("insertions_check")){
  if(nrow(insertions_check) == 0){
    pav_ins_plotting[,9] <- 0
    pav_ins_plotting[,10] <- 0
  }
}

##Create inversion dataframe for plotting (if inversions are present)
if(nrow(inversions) > 0){
  
  make_inv_polygon_forward <- function(start, y_pos, bracket_width){
    tibble(
      x = c(
        as.numeric(start) - size_buffer, 
        as.numeric(start) - size_buffer + bracket_width,
        as.numeric(start) + size_buffer, 
        as.numeric(start) + size_buffer, 
        as.numeric(start) - size_buffer + bracket_width, 
        as.numeric(start) - size_buffer
      ),
      y = c(
        y_pos,
        y_pos, 
        y_pos + 0.02, 
        y_pos + 0.25 - 0.02, 
        y_pos + 0.25, 
        y_pos + 0.25
      )
    )
  }
  
  make_inv_polygon_reverse <- function(end, y_pos, bracket_width){
    tibble(
      x = c(
        as.numeric(end) + size_buffer - bracket_width, 
        as.numeric(end) + size_buffer,
        as.numeric(end) + size_buffer, 
        as.numeric(end) + size_buffer - bracket_width, 
        as.numeric(end) - size_buffer, 
        as.numeric(end) - size_buffer
      ),
      y = c(
        y_pos,
        y_pos, 
        y_pos + 0.25, 
        y_pos + 0.25, 
        y_pos + 0.25 - 0.02, 
        y_pos + 0.02
      )
    )
  }
  
  pav_inv_plotting_forward <- inversions %>%
    dplyr::rowwise() %>%
    dplyr::mutate(polygon = list(make_inv_polygon_forward(start = as.numeric(START), y_pos = -0.15, bracket_width = 0.005 * interval_length))) %>%
    unnest(polygon)
  
  pav_inv_plotting_reverse <- inversions %>%
    dplyr::rowwise() %>%
    dplyr::mutate(polygon = list(make_inv_polygon_reverse(end = as.numeric(END), y_pos = -0.15, bracket_width = 0.005 * interval_length))) %>%
    unnest(polygon)
  
}





##JOSH## plot 3
##First, alter the dataframe to add gene layer information for plotting
assign_layers <- function(N2_df) {
  
  #Make sure to arrange by start coord or this wont work 
  #Initialize column for layer in the dataframe 
  N2_df <- N2_df %>% dplyr::arrange(newstart)
  N2_df$layer <- NA_integer_
  layers <- list()
  
  #For each gene, get the start and stop coords within the interval 
  for (i in  1:nrow(N2_df)){ 
    start <- N2_df$newstart[i]
    end <- N2_df$newend[i]
    
    #Assume the gene has not been placed in a layer 
    placed <- FALSE
    
    #Check each layer to see if the gene overlaps with the last stop coord in that layer
    #If not, place the gene in that layer 
    for (t in seq_along(layers)) {
      last_end <- layers[[t]]
      if (start > last_end + (0.01 * interval_length)) {
        N2_df$layer[i] <- t
        layers[[t]] <- end
        placed <- TRUE
        break
      }
    }
    #If the gene is not placed in any layer, make a new layer, and update the final stop coord in that new layer 
    if (!placed){
      N2_df$layer[i] <- length(layers) + 1
      layers[[length(layers) + 1]] <- end
    }
    
  }
  return(N2_df)
}

#Now, update the dataframe again with information to construct the gene polygon 

make_gene_polygon <- function(start, end, ycenter, strand, width = 0.8, arrow_length){
  if(strand == "+"){
    #Make a dataframe with the x and y coordinates that will make up the gene polygon shape for each row in the main dataframe
    tibble(
      x = c(start, end, end + arrow_length, end, start),
      y = c(ycenter - width/2, ycenter - width/2, ycenter, ycenter + width/2, ycenter + width/2)
    )
  }else{
    #Reverse the polygon so that the arrows points in the other direction (for - strands)
    tibble(
      x = c(start, end, end, start, start - arrow_length),
      y = c(ycenter - width/2, ycenter - width/2, ycenter + width/2, ycenter + width/2, ycenter)
    )
  }
}

#Alter the N2 df with information for layer-wise plotting and gene arrow plotting. Also converts from wide to long format, with data for each polygon point. 
N2_arrowplot <- N2ad %>%
  assign_layers() %>%
  #Flip the order of the layers 
  dplyr::mutate(layer = layer * -1) %>%
  dplyr::mutate(gene_id = row_number()) %>%
  dplyr::rowwise() %>%
  #0.005 * interval length: the maximum distance an arrow can extend past the gene start/stop. Normalizes for the size of the interval selected. 
  dplyr::mutate(polygon = list(make_gene_polygon(newstart, newend, layer, strand, arrow_length = 0.005 * (interval_length)))) %>%
  unnest(polygon) %>%
  #Set y axis back to positive
  dplyr::mutate(y = y + max(abs(layer)) + 1)

#repeat for pseudogenes
N2_pseudogenes_arrowplot <- N2ad_pseudo %>%
  assign_layers() %>%
  #Flip the order of the layers 
  dplyr::mutate(layer = layer * -1) %>%
  dplyr::mutate(gene_id = row_number()) %>%
  dplyr::rowwise() %>%
  #0.005 * interval length: the maximum distance an arrow can extend past the gene start/stop. Normalizes for the size of the interval selected. 
  dplyr::mutate(polygon = list(make_gene_polygon(newstart, newend, layer, strand, arrow_length = 0.005 * (interval_length)))) %>%
  unnest(polygon) %>%
  #Set y axis back to positive
  dplyr::mutate(y = y + max(abs(layer)) + 1)



##Finally, plot the gene models 
plot3 <- ggplot() +
  geom_polygon(data = N2_arrowplot, aes(x = as.numeric(x)/1e6, y = as.numeric(y), group = gene_id, fill = strand, linetype = genetype, alpha = genetype), color = "black") +
  scale_fill_manual(values = c("gray50", "gray90")) + # Watson is gray90 (lightgray)
  #Have to make sure x axes are identical between plot3 and plot4!
  scale_x_continuous(
    limits = c((as.numeric(args_start) - 0.01*(interval_length))/1e6, (as.numeric(args_stop) + 0.01*(interval_length))/1e6)
    ) + 
  scale_alpha_manual(values = c("protein_coding" = 1, "pseudogene" = 0.5)) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),  
    plot.title = element_text(size = 16, face = "bold", hjust = 0.3),
    legend.position = "none",
    panel.background = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )


plot3

##JOSH##


plot4 <- ggplot() +
  geom_rect(data = GATK, aes(xmin = (as.numeric(start) - size_buffer)/1e6, xmax = (as.numeric(end) + size_buffer)/1e6, ymin = 1.05, ymax = 1.3, fill = 'SNVs (GATK)')) + 
  geom_rect(data = deletions, aes(xmin = (as.numeric(START) - size_buffer)/1e6, xmax = (as.numeric(END) + size_buffer)/1e6, ymin = 0.65, ymax = 0.9, fill = 'Deletions')) +
   #geom_point(data = deletions, aes(x = as.numeric(START)/1e6, y = 0.001, fill = 'Deletions')) + 
  geom_rect(data = G_del, aes(xmin = (as.numeric(START) - size_buffer)/1e6, xmax = (as.numeric(END) + size_buffer)/1e6, ymin = 0.65, ymax = 0.9, fill = 'Deletions')) +
   #geom_point(data = G_del, aes(x = as.numeric(START)/1e6, y = 0.001, fill = 'Deletions')) + 
  geom_polygon(data = gatk_ins_plotting, aes(x = as.numeric(x)/1e6, y = as.numeric(y), group = POS, fill = 'Insertions')) + #Dont need to add size buffer for insertions, already did that earlier 
  geom_polygon(data = pav_ins_plotting, aes(x = as.numeric(x)/1e6, y = as.numeric(y), group = POS, fill = 'Insertions')) + 
  geom_text(data = pav_ins_plotting, aes(x = as.numeric(START)/1e6, y = 0.5 + 0.1, label = SV_length), size = 1.5) + #Label large PAV insertions
  geom_point(data = data.frame(x = c(16187171/1e6,16190239/1e6), y = c(0,0)), mapping = aes(x = x, y = y)) + 
  
  #Making sure x axes are consistent! 
  scale_x_continuous(
    name = "Genomic position (Mb)", 
    labels = scales::number_format(scale = 1, accuracy = 0.01),
    limits = c((as.numeric(args_start) - 0.01*(interval_length))/1e6, (as.numeric(args_stop) + 0.01*(interval_length))/1e6)) +
  scale_fill_manual(
    name = "",
    values = c("SNVs (GATK)" = "black", "Insertions" = "blue", "Deletions" = "red"),
    breaks = c("SNPs (GATK)", "Insertions", "Deletions")
  ) + 
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 10, color = "black"),  
    axis.title.x = element_text(size = 12, color = "black"),
    panel.grid = element_blank(), 
    plot.background = element_rect(fill = "white", color = NA),  
    axis.line.x = element_line(),
    legend.position = "none",
    panel.background = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

plot4

##Plot SNP annotations on a separate plot 
if(nrow(bcsq_annot_df) > 0){
  library(viridis)
  
  plot4_annot <- ggplot() +
    geom_rect(data = GATK_annot, aes(xmin = (as.numeric(start) - size_buffer)/1e6, xmax = (as.numeric(end) + size_buffer)/1e6, ymin = 1.05, ymax = 1.3, fill = bcsq_annot)) + 
    # geom_rect(data = deletions, aes(xmin = (as.numeric(START) - size_buffer)/1e6, xmax = (as.numeric(END) + size_buffer)/1e6, ymin = 0.65, ymax = 0.9), fill = 'red') +
    # #geom_point(data = deletions, aes(x = as.numeric(START)/1e6, y = 0.001, fill = 'Deletions')) + 
    # geom_rect(data = G_del, aes(xmin = (as.numeric(START) - size_buffer)/1e6, xmax = (as.numeric(END) + size_buffer)/1e6, ymin = 0.65, ymax = 0.9), fill = 'red') +
    # #geom_point(data = G_del, aes(x = as.numeric(START)/1e6, y = 0.001, fill = 'Deletions')) + 
    # geom_polygon(data = gatk_ins_plotting, aes(x = as.numeric(x)/1e6, y = as.numeric(y), group = POS), fill = 'blue') + #Dont need to add size buffer for insertions, already did that earlier 
    # geom_polygon(data = pav_ins_plotting, aes(x = as.numeric(x)/1e6, y = as.numeric(y), group = POS), fill = 'blue') + 
    # geom_text(data = pav_ins_plotting, aes(x = as.numeric(START)/1e6, y = 0.5 + 0.1, label = SV_length), size = 1.5) + #Label large PAV insertions
     scale_fill_manual(
      values = c(
        "stop gained"    = "black",  
        "missense"       = "#D62728",  
        "3 prime UTR"    = "#FFD92F",  
        "5 prime UTR"    = "#E6AC00",  
        "noncoding"      = "grey50",  
        "synonymous"     = "#17BECF",  
        "intron"         = "#1F77B4",  
        "splice region"  = "#66C2A5"
      ), na.value = "grey90"
    ) +
  
  
  
  #Making sure x axes are consistent! 
    scale_x_continuous(
      name = "Genomic position (Mb)", 
      labels = scales::number_format(scale = 1, accuracy = 0.01),
      limits = c((as.numeric(args_start) - 0.01*(interval_length))/1e6, (as.numeric(args_stop) + 0.01*(interval_length))/1e6)) +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 10, color = "black"),  
      axis.title.x = element_text(size = 12, color = "black"),
      panel.grid = element_blank(), 
      plot.background = element_rect(fill = "white", color = NA),  
      axis.line.x = element_line(),
     legend.position = "bottom",
      panel.background = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    ) +
    labs(fill = NULL) 

  plot4_annot
}




#Make a slightly different plot if there are inversions

if(nrow(inversions) > 0){

  plot4 <- ggplot() +
    geom_rect(data = GATK, aes(xmin = (as.numeric(start) - size_buffer)/1e6, xmax = (as.numeric(end) + size_buffer)/1e6, ymin = 1.05, ymax = 1.3, fill = 'SNVs (GATK)')) + 
    geom_rect(data = deletions, aes(xmin = (as.numeric(START) - size_buffer)/1e6, xmax = (as.numeric(END) + size_buffer)/1e6, ymin = 0.65, ymax = 0.9, fill = 'Deletions')) +
    geom_rect(data = G_del, aes(xmin = (as.numeric(START) - size_buffer)/1e6, xmax = (as.numeric(END) + size_buffer)/1e6, ymin = 0.65, ymax = 0.9, fill = 'Deletions')) +
    geom_polygon(data = gatk_ins_plotting, aes(x = as.numeric(x)/1e6, y = as.numeric(y), group = POS, fill = 'Insertions')) + #Dont need to add size buffer for insertions, already did that earlier 
    geom_polygon(data = pav_ins_plotting, aes(x = as.numeric(x)/1e6, y = as.numeric(y), group = POS, fill = 'Insertions')) + 
    geom_text(data = pav_ins_plotting, aes(x = as.numeric(START)/1e6, y = 0.5 + 0.1, label = SV_length), size = 1.5) + #Label large PAV insertions
    #Add inversions to the plot    
    geom_polygon(data = pav_inv_plotting_forward, aes(x = as.numeric(x)/1e6, y = as.numeric(y), group = POS, fill = 'black'), color = "black") + 
    geom_polygon(data = pav_inv_plotting_reverse, aes(x = as.numeric(x)/1e6, y = as.numeric(y), group = POS, fill = 'black'), color = "black") + 
    scale_x_continuous(
      name = "Genomic position (Mb)", 
      labels = scales::number_format(scale = 1, accuracy = 0.01),
      limits = c((as.numeric(args_start) - 0.01*(interval_length))/1e6, (as.numeric(args_stop) + 0.01*(interval_length))/1e6)) +
    scale_fill_manual(
      name = "",
      values = c("SNVs (GATK)" = "black", "Insertions" = "blue", "Deletions" = "red"),
      breaks = c("SNPs (GATK)", "Insertions", "Deletions")
    ) + 
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 10, color = "black"),  
      axis.title.x = element_text(size = 12, color = "black"),
      panel.grid = element_blank(), 
      plot.background = element_rect(fill = "white", color = NA),  
      axis.line.x = element_line(),
      legend.position = "none",
      panel.background = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )
  
}





combined_plot2 <- plot_grid(plot3, plot4, ncol = 1, align = "v", rel_heights = c(0.35, 0.65))

combined_plot2

if(nrow(bcsq_annot_df) > 0) {
  combined_plot2_annot <- plot_grid(plot3, plot4_annot, ncol = 1, align = "v", rel_heights = c(0.35, 0.65))

  combined_plot2_annot
}
  
if(nrow(inversions) > 0){
  combined_plot2 <- plot_grid(plot3, plot4, ncol = 1, align = "v", rel_heights = c(0.25, 0.75))
  
}  

plot_file <- paste("/vast/eande106/projects/Josh/SV_calling/processed_data/vcf_filtering", 
                    args_species, paste(
                      "variant_plot",
                      "_",
                      args_species, 
                      "_",
                      args_chrom, 
                      "_",
                      args_start, 
                      "_",
                      args_stop, 
                      "_",
                      args_strain,
                      ".jpg",
                      sep = ""
                    ), sep = "/"
)

plot_file_annotated <- paste("/vast/eande106/projects/Josh/SV_calling/processed_data/vcf_filtering", 
                   args_species, paste(
                     "variant_plot_annotated",
                     "_",
                     args_species, 
                     "_",
                     args_chrom, 
                     "_",
                     args_start, 
                     "_",
                     args_stop, 
                     "_",
                     args_strain,
                     ".jpg",
                     sep = ""
                   ), sep = "/"
)



ggsave(plot_file, combined_plot2, dpi=600, width = 7.5, height = 1.8)

if(nrow(bcsq_annot_df) > 0){
  ggsave(plot_file_annotated, plot4_annot, dpi=600, width = 7.5, height = 3) #1.8)
}

##make df report for deletions (including overlapping genes)
deletions_sv_report <- deletions %>%
  dplyr::rename("Variant" = SV) %>%
  dplyr::rename("Var_length" = SV_length)
deletions_sv_report$called_by <- "PAV"

deletions_gatk_report <- G_del
deletions_gatk_report$called_by <- "GATK"

deletions_report <- deletions_sv_report %>%
  rbind(deletions_gatk_report) %>%
  dplyr::arrange(POS) %>%
  fuzzyjoin::interval_left_join(N2ad, by = c("START" = "newstart", "END" = "newend")) %>%
  dplyr::rename("gene" = N2) %>%
  dplyr::select(c("CHROM", "Var_length", "START", "END", "called_by", "gene", "REF", "ALT"))

##same for insertions
insertions_sv_report <- insertions %>%
  dplyr::rename("Variant" = SV) %>%
  dplyr::rename("Var_length" = SV_length)
insertions_sv_report$called_by <- "PAV"

insertions_gatk_report <- G_ins
insertions_gatk_report$called_by <- "GATK"

insertions_report <- insertions_sv_report %>%
  rbind(insertions_gatk_report) %>%
  dplyr::mutate("END" = START + 1) %>%
  dplyr::arrange(POS) %>%
  fuzzyjoin::interval_left_join(N2ad, by = c("START" = "newstart", "END" = "newend")) %>%
  dplyr::rename("gene" = N2) %>%
  dplyr::select(c("CHROM", "POS", "Var_length", "called_by", "gene", "REF", "ALT"))

#Make SNP report
SNP_report <- GATK_annot %>%
  fuzzyjoin::interval_left_join(N2ad, by = c("start" = "newstart", "end" = "newend")) %>%
  dplyr::rename("Overlapping gene" = N2) %>%
  dplyr::select("CHROM", "start", "REF", "ALT", "Overlapping gene", "bcsq_annot", "annovar_annot", "vep_annot") %>%
  dplyr::rename("POS" = start) %>%
  dplyr::rename("BCSQ Annotation" = bcsq_annot) %>%
  dplyr::rename("ANNOVAR Annotation" = annovar_annot) %>%
  dplyr::rename("VEP Annotation" = vep_annot)

#save.image(file = "/vast/eande106/projects/Josh/SV_calling/processed_data/variantplotting.RData")


