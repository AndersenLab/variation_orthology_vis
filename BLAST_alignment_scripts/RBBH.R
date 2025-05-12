library(ape)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)


##Create gene:transcript table for N2
N2gff <- ape::read.gff("/vast/eande106/data/c_elegans/genomes/PRJNA13758/WS283/csq/c_elegans.PRJNA13758.WS283.csq.gff3")
N2_txpt <- N2gff %>%
  dplyr::filter(type == "mRNA") %>%
  dplyr::select(attributes) %>%
  tidyr::separate(attributes, into = as.character(c(1:16)), sep = "[:;=]") %>%
  dplyr::select(c("3", "6"))
 
colnames(N2_txpt) <- c("transcript", "N2gene")

#Open wi->N2 protein blast
wi_protb <- readr::read_tsv(
  "/vast/eande106/projects/Josh/alignment_BLAST/blastp/CB4856.braker.protein_wi_blast_N2lib.pb.txt",
  col_names = as.character(c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                             "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
  )


wi_protb_split <- wi_protb %>%
  #remove nonsignficant hits 
  filter(evalue < 1e-5) %>%
  #Split query transcript column into gene and transcript 
  tidyr::separate(qseqid, into = c("qgene", "qtranscript"), sep = "[.]") %>%
#Split subject transcript column into gene and transcript
  tidyr::separate(sseqid, into = c("txpt", "stranscript"), sep = "[:]") %>%
  select(!txpt)


gene_groups <- wi_protb_split %>%
  #Add N2 gene names via left-join
  left_join(N2_txpt, by = c("stranscript" = "transcript")) %>%
  #Group by query gene
  group_by(qgene) %>%
#Find the best transcript hit for each wi gene 
  arrange(evalue,  desc(bitscore), .by_group = TRUE) %>%
  slice_head(n = 1)


#Read in the N2-> wi blastp
N2_protb <- readr::read_tsv(
  "/vast/eande106/projects/Josh/alignment_BLAST/blastp/CB4856.braker.protein_N2_blast_wilib.pb.txt",
  col_names = as.character(c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                             "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
)

N2_protb_split <- N2_protb %>%
  #Remove nonsignificant hits 
  filter(evalue < 1e-5) %>%
  #Split query transcript column into gene and transcript 
  tidyr::separate(qseqid, into = c("txpt", "qtranscript"), sep = "[:]") %>%
  select(!txpt) %>%
#Split subject transcript column into gene and transcript
  tidyr::separate(sseqid, into = c("sgene", "stranscript"), sep = "[.]") 


N2_blast <- N2_protb_split %>%
  #left join by N2 transcript name to include N2 genes
  left_join(N2_txpt, by = c("qtranscript" = "transcript")) %>%
  group_by(N2gene) %>%
 #Find the best blast N2 transcript hit for each gene 
  arrange(evalue,  desc(bitscore), .by_group = TRUE) %>%
  slice_head(n = 1)

#Finally, perform RBBH
RBBH <- gene_groups %>%
  left_join(N2_blast, by = c("N2gene" = "N2gene")) %>%
  filter(qgene == sgene)
