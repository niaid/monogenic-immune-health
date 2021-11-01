library(tidyverse)

#setwd("/hpcdata/sg/sg_data/PROJECTS/Monogenic_Project")

JIVE.PC.IN.PATH <- snakemake@input[["jive_pca"]]#"Integration_output/jive/subject/prcomp_list.rds"
JIVE.IN.PATH <- snakemake@input[["jive"]]#"Integration_output/jive/subject/jive.rds"

TABLE.OUT.PATH <- snakemake@output[[1]]#"Paper_1_Figures/Figure_3_Tables/jive_pcs.csv"


prcomp.list <- readRDS(JIVE.PC.IN.PATH)
jive <- readRDS(JIVE.IN.PATH)
pdat <- jive$pdat

pdat <- pdat %>% select(patient_id, condition)

joint <- prcomp.list$joint$x[, 1:2]
colnames(joint) <- paste0("j", colnames(joint))

gene <- prcomp.list$array.ind$x[, 1:2]
colnames(gene) <- paste0("transcriptome_i", colnames(gene))

protein <- prcomp.list$soma.ind$x[, 1:2]
colnames(protein) <- paste0("proteome_i", colnames(protein))

pcs <- do.call(cbind, list(joint, gene, protein))


stopifnot(identical(rownames(pcs), pdat$patient_id))

dat <- 
  bind_cols(pdat, as.data.frame(pcs))


write_csv(dat, TABLE.OUT.PATH)
