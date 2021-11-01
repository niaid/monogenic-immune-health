
# Load relevant libraries
library(Biobase)

# Source utilities

if(!exists("snakemake")){
  setwd("../../..")
  source("scripts/util/paper/parse_snakemake.R")
  parse_snakemake(rule = "methods_stats_patients_missing_cbc")
}

# Set globals
## the tbnk training eset (prior to cleaning)
TBNKS.IN.PATH = snakemake@input[[1]]#'Data/TBNK/processed/tbnk_eset_training.rds'

## the cleaned up sample-level tbnk training eset
TBNKS.SAMPLE.LEVEL.OUT.PATH = snakemake@output[[1]]#'Data/TBNK/data_analysis_ready/tbnk_training_sample_level_eset.rds'
## the cleaned up subject-level tbnk training eset
TBNKS.SUBJECT.LEVEL.OUT.PATH = snakemake@output[[2]]#'Data/TBNK/data_analysis_ready/tbnk_training_subject_level_eset.rds'

# Load data
tbnks = readRDS(TBNKS.IN.PATH)

# Remove unwanted samples for the following reasons
# V316 & V318: Lymphocyte subpopulations reported in TBNKs do not sum to total lymphocytes
# V40: Lymphocyte percentages reported in TBNKs are not equal to absolute population counts divided by total lymphocytes
outliers = c('V40', 'V316','V318')
tbnks = tbnks[, ! tbnks$visit_id %in% outliers]

mat <- exprs(tbnks)
sum_na <- rowSums(is.na(mat))

sum_na

#mat["nucleated_rbc_abs", ]
