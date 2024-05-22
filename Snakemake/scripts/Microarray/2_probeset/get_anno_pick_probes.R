# Picks the probeset most correlated with PC1 of the probesets for all of the probes. 
# Uses utility function from old microarray pipeline
# The output is the PICKED.PROBES.OUT.PATH file where the picked probes are saved

library(plyr) #need for the pick probeset function
library(dplyr)
library(tidyr)
library(GEOquery)

#sources the function that will find probeset most correlated with PC1 of all probesets for each gene
source("scripts/util/Processing/pick_probeset.R")

# Config ------------------------------------------------------
#input files
#path to microarray data
ARRAY.IN.PATH <- snakemake@input[["eset"]]#"Data/Microarray/raw/eset_rma_with_pData.rds"
ANNOTATION.IN.PATH <- snakemake@input[["probe_anno"]]#"Data/Microarray/probeset/pre_downloaded_ann/probe_annotations_full.csv"

#outfiles
PROBEMAP.OUT.PATH <- snakemake@output[["probemap"]]#"Data/Microarray/probeset/output/probe_annotations.txt" #saves probemap file that is input to the pick probeset function
PICKED.PROBES.OUT.PATH <- snakemake@output[["picked_probes"]]#"Data/Microarray/probeset/output/picked_probes.txt"

#if(!dir.exists(dirname(PICKED.PROBES.OUT.PATH))) dir.create(dirname(PICKED.PROBES.OUT.PATH))

#load data ------------------------------------------------------
eset <- readRDS(ARRAY.IN.PATH) 
annotation <- read.csv(ANNOTATION.IN.PATH, header = TRUE)

#write the probemap file that maps genes to probes
out <- 
  annotation %>% 
  select(ID, Gene.symbol) %>% # keeps only ID and gene columns
  rename(gene = Gene.symbol) %>% # rename column
  mutate(ID = as.character(ID), gene = as.character(gene)) %>%  # Make sure they are character vectors
  filter(!grepl("///",gene) & gene != "")              # Filter out genes that map to more than one ID

write.table(out, file = PROBEMAP.OUT.PATH, sep="\t", quote = F, col.names = T, row.names = F)

# Call pick probeset function. Will write the picked probes to file

pick.probeset(eset, PROBEMAP.OUT.PATH, PICKED.PROBES.OUT.PATH) # generates file.map.pc1
