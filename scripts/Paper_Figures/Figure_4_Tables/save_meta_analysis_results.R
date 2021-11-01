# Load library
library(tidyverse)

# Set paths
#setwd("/hpcdata/sg/sg_data/PROJECTS/Monogenic_Project")
VALIDATION.RESULTS.IN.PATH = snakemake@input[["overall_res"]]#"Reference/jamboree/analysis_output/results/jamboree_z_score_results.RDS"
POOLED.TABLE.OUT.PATH = snakemake@output[["pooled_tab"]]#'Paper_1_Figures/Figure_4_Tables/figure_4_meta_analysis_table.txt'
EFFSIZE.TABLE.OUT.PATH = snakemake@output[["effect_sizes"]]#'Paper_1_Figures/Figure_4_Tables/figure_4_meta_analysis_effsize_table.txt'

source("scripts/util/Plotting/tbnk_featurename_replace.R")
# Load data
scores = readRDS(VALIDATION.RESULTS.IN.PATH)

rownames(scores$metaAnalysis$pooledResults)

pooled <- scores$pooledResults %>%
        rownames_to_column(var = "feature")

pooled <- pooled %>% filter(feature != "microarray.classifier") %>%
        mutate(feature = gsub("PC1", "jPC1", feature)) %>%
        mutate(feature = gsub("somalogic\\.grey", "serum", feature)) %>%
        mutate(feature = gsub("somalogic\\.modules\\.purple", "PM2", feature)) %>%
        mutate(feature = gsub("tbnks\\.", "", feature)) %>%
        mutate(feature = gsub("healthy\\.index", "Immune Health Metric", feature)) %>%
        mutate(feature = replace_tbnk_names(feature))

to_dat <- function(mat, value_col_name){
  dat <- as.data.frame(mat) %>% rownames_to_column(var = "feature") %>%
          gather(key = "study", value = "value", -feature)

  colnames(dat)[colnames(dat)== "value"] <- value_col_name
  dat
}

effect_sizes <- left_join(
  to_dat(scores$datasetEffectSizes, "effectSize"),
  to_dat(scores$datasetEffectSizes, "effectSizeStandardError")
)

effect_sizes <- effect_sizes %>%
        mutate(feature = gsub("PC1", "jPC1", feature)) %>%
        mutate(feature = gsub("somalogic\\.grey", "serum", feature)) %>%
        mutate(feature = gsub("somalogic\\.modules\\.purple", "PM2", feature)) %>%
        mutate(feature = gsub("tbnks\\.", "", feature)) %>%
        mutate(feature = gsub("healthy\\.index", "Immune Health Metric", feature)) %>%
        mutate(feature = replace_tbnk_names(feature))



write_csv(pooled, POOLED.TABLE.OUT.PATH)
write_csv(effect_sizes, EFFSIZE.TABLE.OUT.PATH)
