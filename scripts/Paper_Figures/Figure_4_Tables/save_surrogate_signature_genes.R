library(tidyverse)

#setwd("/hpcdata/sg/sg_data/PROJECTS/Monogenic_Project")

SIGS.IN.PATH <- snakemake@input[[1]]#"Classification/transcriptional_surrogates/surrogate_signatures.RDS"

TABLE.OUT.PATH <- snakemake@output[[1]]#"Paper_1_Figures/Figure_4_Tables/surrogate_sig_genes.csv"

source("scripts/util/Plotting/tbnk_featurename_replace.R")

sigs <- readRDS(SIGS.IN.PATH)

dat_list <- lapply(sigs, function(x){
  sig_members <- sapply(x, paste, collapse = " ")
  data.frame(direction = names(x), gene_symbols = sig_members)
})

dat <- bind_rows(dat_list, .id = "signature")


dat <- dat %>% filter(signature != "microarray.classifier") %>%
        mutate(signature = gsub("PC1", "jPC1", signature)) %>%
        mutate(signature = gsub("somalogic\\.grey", "serum", signature)) %>%
        mutate(signature = gsub("somalogic\\.modules\\.purple", "PM2", signature)) %>%
        mutate(signature = gsub("tbnks\\.", "", signature)) %>%
        mutate(signature = gsub("healthy\\.index", "Immune Health Metric", signature))

dat$signature <- replace_tbnk_names(dat$signature)


write_csv(dat, TABLE.OUT.PATH)
