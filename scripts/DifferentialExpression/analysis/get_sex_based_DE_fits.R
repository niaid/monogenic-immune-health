## Fits differential expression linear models using limma. Model is feature ~ age + sex*condition + visit_type.
## Visit type corresponds to whether the patient was known to be sick upon visiting.

# Load libraries
library(limma)
library(Biobase)
library(dplyr)
source('scripts/util/DifferentialExpression/limma.R')

# Set seed		
set.seed(19020)

# Set globals
## Path to esets used for fitting the model
ESET.IN.PATHS = list(
  somalogic.modules = snakemake@input[[1]],#'Data/Somalogic/analysis_output/wgcna_results/scores_sample_level.rds',
  somalogic.features = snakemake@input[[2]],#'Data/Somalogic/data_analysis_ready/analysis_ready_sample_level_training_somalogic.rds',
  microarray.modules = snakemake@input[[3]],#'Data/Microarray/analysis_output/WGCNA/array_sample_scores.rds',
  microarray.features = snakemake@input[[4]],#'Data/Microarray/data_analysis_ready/eset_batch_training_sample.rds',
  tbnks.features = snakemake@input[[5]]#'Data/TBNK/data_analysis_ready/tbnk_training_sample_level_eset.rds'
)

## Where to save the fitted sex-linked DE 
limma models
FIT.OUT.PATHS = list(
  somalogic.modules = snakemake@output[[1]],#'Data/Somalogic/analysis_output/sex_related_de_signatures/somalogic_modules_DE_sex_linked_fit.rds',
  somalogic.features = snakemake@output[[2]],#'Data/Somalogic/analysis_output/sex_related_de_signatures/somalogic_features_DE_sex_linked_fit.rds',
  microarray.modules = snakemake@output[[3]],#'Data/Microarray/analysis_output/sex_related_de_signatures/microarray_modules_DE_sex_linked_fit.rds',
  microarray.features = snakemake@output[[4]],#'Data/Microarray/analysis_output/sex_related_de_signatures/microarray_features_DE_sex_linked_fit.rds',
  tbnks.features = snakemake@output[[5]]#'Data/TBNK/analysis_output/sex_related_de_signatures/tbnks_features_DE_sex_linked_fit.rds'
)

# Load esets
esets = lapply(ESET.IN.PATHS, readRDS)

# Instantiate a function to get stats associated with stable features of an eset
get_fit = function(eset, condition) {
  
  eset = eset[, eset$condition %in% c(condition, 'Healthy')]
  
  # Scale eset
  exprs(eset) = t(scale(t(exprs(eset))))

  # Get design matrix
  design = make_sex_linked_design(eset)
  
  # Fit limma model
  fit = fit_limma(eset, design)
  
  return(fit)
}

get_fits = function(eset) {
  conditions = c('STAT1 GOF', '47CGD', 'Job')
  fits = lapply(conditions, get_fit, eset = eset)
  names(fits) = conditions
  return(fits)
}

# Apply the function over all esets
fits = lapply(esets, get_fits)

# Save fits
mapply(function(fit, out.path) {
  saveRDS(fit, out.path)
}, fits, FIT.OUT.PATHS)