library(plyr) # ALWAYS LOAD BEFORE DPLYR
library(dplyr)
library(ggplot2)
library(readr)
library(IRanges)
library(valr)
library(data.table)




# Change VCF to PAV INDEL calls
#This file comes from pavmerge_trim.sh
vcf <- readr::read_tsv("/vast/eande106/projects/Josh/SV_calling/raw_data/c_elegans/pav_vcfs/pav_indels_merged.tsv") 

#Filter out JU2617
vfilt <- vcf %>%
  dplyr::select(!JU2617) %>%
  dplyr::filter(`#CHROM` != "MtDNA") %>% 
  dplyr::rename(CHROM=`#CHROM`)



pos_bed <- vfilt %>% 
  dplyr::select(CHROM, POS, SV) %>%
  dplyr::rename(start=POS) %>%
  dplyr::mutate(end=start+1) %>% 
  dplyr::rename(chrom=CHROM) %>%
  dplyr::select(chrom, start, end, SV)



ggplot() + 
  # geom_rect(data = all_collapsed %>% dplyr::rename(CHROM = chrom), aes(xmin = start/1e6, xmax = end/1e6, ymin = 0, ymax = 1.75, fill = "hdr")) + 
  geom_rect(data= vfilt %>% dplyr::filter(SV == "DEL"), aes(xmin = POS/1e6, xmax = POS/1e6 + 0.001, ymin = 1, ymax = 1.5, fill = "DEL")) +
  geom_rect(data = vfilt %>% dplyr::filter(SV == "INS"), aes(xmin = POS/1e6, xmax = POS/1e6 + 0.001, ymin = 0.25, ymax = 0.75, fill = "INS")) +
  # geom_rect(data = )
  facet_wrap(~CHROM, nrow=6, scales = "free_x") + 
  scale_fill_manual(
    name = "",
    values = c("DEL" = "red", 'INS' = 'blue', 'hdr' = 'gray70'),
    breaks = c("Deletions", "Insertion")
  ) + 
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 10, color = "black"),  
    axis.title.x = element_text("Genome position (Mb)", size = 12, color = "black"),
    panel.grid = element_blank(), 
    plot.background = element_rect(fill = "white", color = NA),  
    axis.line.x = element_line(),
    legend.position = "none",
    strip.text = element_text(size = 14, color = 'black'),
    panel.background = element_blank(),
    # plot.margin = margin(0, 0, 0, 0)
  )



vfilt_longer <- vfilt %>%
  dplyr::select(-REF,-ALT, -SV_length) %>%
  pivot_longer(
    cols = -c(CHROM, POS, SV),  # Keep chrom, pos, and info as identifiers
    names_to = "sample",          # New column for sample names
    values_to = "genotype"           # New column for sample values
    ) %>%
  dplyr::mutate(genotype=ifelse(genotype=="./.",0,1)) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(n_alt=sum(genotype)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(n_alt)) %>%
  dplyr::mutate(rid=rleid(n_alt)) %>%
  dplyr::filter(genotype > 0)
#save as text file 
write.table(vfilt_longer, file = "/vast/eande106/projects/Josh/SV_calling/raw_data/c_elegans/vfilt_longer.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


### Then you will need to intersect with HDRs for a particular strain and add a column to indicate if that variant is found in a HDR for that particular strain
INDEL_plt <- ggplot(vfilt_longer) + geom_rect(aes(xmin=(POS-500)/1e6,xmax=(POS+500)/1e6,ymin=rid+0.7,ymax=rid-0.7,fill=SV)) + 
  facet_wrap(~CHROM,nrow=1,scales = 'free_x') +
  scale_fill_manual(values = c("DEL" = "red", "INS" = "blue")) +
  theme(
    panel.background = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.border = element_rect(fill = NA),
    axis.title = element_text(size = 13),
    legend.position = "none") +
  xlab("Genome position (Mb)") +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) + 
  ylab("114 WI genomes")
INDEL_plt
ggsave("/vast/eande106/projects/Josh/SV_calling/processed_data/pav_indel_plot.jpg", INDEL_plt, dpi = 300, width = 14, height = 10)


#Get HDR dataframe
hdrs <- readr::read_tsv("/vast/eande106/projects/Josh/SV_calling/raw_data/c_elegans/114HDRs.tsv", col_names = FALSE)
colnames(hdrs) <- c("chrom", "start", "stop", "strain")
#save as text file 
write.table(hdrs, file = "/vast/eande106/projects/Josh/SV_calling/raw_data/c_elegans/hdrs.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##Run the pav_hdr_intersect.sh script (uses bedtools intersect)##

#This gives non-hdr variants
vfilt_longer_nohdrs <- read_tsv("/vast/eande106/projects/Josh/SV_calling/raw_data/c_elegans/vfilt_longer_no_hdrs.bed", col_names = FALSE)
vfilt_longer_nohdrs <- vfilt_longer_nohdrs %>%
  tidyr::separate(X1, into = c("chrom", "strain"), sep = "[_]")

non_hdr_rows <- vfilt_longer_nohdrs$X4

#pull non-hdr variants from the original filtered dataframe
vfilt_longer_nohdrs_plotting <- vfilt_longer[non_hdr_rows, ]


# Plot only INDELs not in a HDR for each strain - code will look something like this
INDEL_plt_noHDR <- ggplot(vfilt_longer_nohdrs_plotting) + geom_rect(aes(xmin=(POS-500)/1e6,xmax=(POS+500)/1e6,ymin=rid+0.7,ymax=rid-0.7,fill=SV)) + 
  facet_wrap(~CHROM,nrow=1,scales = 'free_x') +
  scale_fill_manual(values = c("DEL" = "red", "INS" = "blue")) +
  theme(
    panel.background = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.border = element_rect(fill = NA),
    axis.title = element_text(size = 13),
    legend.position = "none") +
  xlab("Genome position (Mb)") +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) + 
  ylab("114 WI genomes")
INDEL_plt_noHDR

ggsave("/vast/eande106/projects/Josh/SV_calling/processed_data/pav_indel_plot_nohdrs.jpg", INDEL_plt_noHDR, dpi = 300, width = 14, height = 10)

