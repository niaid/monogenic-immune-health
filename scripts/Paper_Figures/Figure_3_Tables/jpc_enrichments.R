library(tidyverse)

#setwd("/hpcdata/sg/sg_data/PROJECTS/Monogenic_Project")

ENRICH.IN.PATH <- snakemake@input[["enrich"]]#"Integration_output/jive/subject/pc_enrich_dat_camera.rds"

TABLE.OUT.PATH <- snakemake@output[[1]]#"Paper_1_Figures/Figure_3_Tables/jive_pc_enrichment.csv"

all.dat <- readRDS(ENRICH.IN.PATH)
all.dat <- all.dat %>% filter(geneset.db != "tiss.general")

all.dat <- all.dat %>%
        filter(PValue < .05) %>%
        mutate(pca.data = replace(pca.data, pca.data == "array.ind", "transcriptome_individual")) %>%
        mutate(pca.data = replace(pca.data, pca.data == "soma.ind", "proteome_individual")) %>%
        filter(PC %in% c("PC1", "PC2"))

all.dat$PC[all.dat$pca.data == "joint"] <- 
        paste0("j", all.dat$PC[all.dat$pca.data == "joint"])

all.dat$PC[all.dat$pca.data == "transcriptome_individual"] <- 
        paste0("transcriptome_i", all.dat$PC[all.dat$pca.data == "transcriptome_individual"])

all.dat$PC[all.dat$pca.data == "proteome_individual"] <- 
        paste0("proteome_i", all.dat$PC[all.dat$pca.data == "proteome_individual"])


all.dat <- all.dat %>%
        mutate(in.data = replace(in.data, in.data == "array", "transcriptome")) %>%
        mutate(in.data = replace(in.data, in.data == "soma", "proteome"))

all.dat <- all.dat %>%
        rename(enrichment.input.data = in.data) %>%
        select(-pca.data) %>%
        select(geneset, everything())

write_csv(all.dat, TABLE.OUT.PATH)
