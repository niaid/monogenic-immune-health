library(tidyverse)

#setwd("/hpcdata/sg/sg_data/PROJECTS/Monogenic_Project")

HI.IN.PATH = snakemake@input[["healthy_index"]]#'Classification/results/healthy_rf_results_all.RDS'
META.IN.PATH = snakemake@input[["meta"]]#'Classification/healthy_random_forest_sample_meta_data_all.RDS'

TABLE.OUT.PATH <- snakemake@output[[1]]#"Paper_1_Figures/Figure_4_Tables/hi_results_full_mod.csv"

results = readRDS(HI.IN.PATH)
meta = readRDS(META.IN.PATH)


## Create a data frame with the HI and condition for each subject
df = data.frame(
                patient_id = rownames(results),
                healthy.index = results$all.modules.plus.grey.with.tbnks, 
                condition = meta$condition)


write_csv(df, TABLE.OUT.PATH)
