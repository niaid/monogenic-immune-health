library(tidyverse)
library(Biobase)
library(BiocGenerics)

#setwd("/hpcdata/sg/sg_data/PROJECTS/Monogenic_Project")
SIG.IN.PATH <- snakemake@input[["sig"]]#"Classification/proteomic_surrogates/healthy.index.plasma.surrogate.somaId.RDS"
ESET.IN.PATH <- snakemake@input[["eset"]]#"Data/Somalogic/data_analysis_ready/analysis_ready_sample_level_training_somalogic.rds"

TABLE.OUT.PATH <- snakemake@output[[1]]#"Paper_1_Figures/Figure_5_Tables/proteomic_surrogate_ihm.csv"

sig <- readRDS(SIG.IN.PATH)
eset <- readRDS(ESET.IN.PATH)


featdata <- featureData(eset)@data %>%
        rename(feature = Target) 


sig_dat <- lapply(sig, function(x){
  data.frame(SomaId= x, stringsAsFactors = FALSE)
}) %>% bind_rows(.id = "direction")

sig_dat <- sig_dat %>% left_join(featdata)

sig_dat <- sig_dat %>% 
        select(-Dilution) %>%
        select(feature, everything())

write_csv(sig_dat, TABLE.OUT.PATH)
