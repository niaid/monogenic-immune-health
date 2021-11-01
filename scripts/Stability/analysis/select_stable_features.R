## Subset relevant sample and subject level esets to the stable features from the simple variance partition model with
## just patient and residual covariates (no medication or condition)

# Load variance partition library
library(variancePartition)
library(Biobase)

# We set the global paths
## Sample level esets
ESET.SAMPLE.LEVEL.IN.PATHS = list(
  somalogic.features = snakemake@input[[1]],#'Data/Somalogic/data_analysis_ready/analysis_ready_sample_level_training_somalogic.rds',
  somalogic.modules = snakemake@input[[2]],#'Data/Somalogic/analysis_output/wgcna_results/scores_sample_level.rds',
  microarray.features = snakemake@input[[3]],#'Data/Microarray/data_analysis_ready/eset_batch_training_sample.rds',
  microarray.modules = snakemake@input[[4]],#'Data/Microarray/analysis_output/WGCNA/array_sample_scores.rds',
  tbnks = snakemake@input[[5]]#'Data/TBNK/data_analysis_ready/tbnk_training_sample_level_eset.rds'
)

## Subject level esets
ESET.SUBJECT.LEVEL.IN.PATHS = list(
  somalogic.features = snakemake@input[[6]],#'Data/Somalogic/data_analysis_ready/analysis_ready_subject_level_training_somalogic.rds',
  somalogic.modules = snakemake@input[[7]],#'Data/Somalogic/analysis_output/wgcna_results/scores_subject_level.rds',
  microarray.features = snakemake@input[[8]],#Data/Microarray/data_analysis_ready/eset_batch_training_subject.rds',
  microarray.modules = snakemake@input[[9]],#'Data/Microarray/analysis_output/WGCNA/array_subject_scores.rds',
  tbnks = snakemake@input[[10]]#'Data/TBNK/data_analysis_ready/tbnk_training_subject_level_eset.rds'
)

## Variance partitions into patient-explained variance and residual variance
VP.IN.PATHS = list(
  somalogic.features = snakemake@input[[11]],#'Data/Somalogic/analysis_output/stability/somalogic_features_standard_vp.rds',
  somalogic.modules = snakemake@input[[12]],#'Data/Somalogic/analysis_output/stability/somalogic_modules_standard_vp.rds',
  microarray.features = snakemake@input[[13]],#'Data/Microarray/analysis_output/stability/microarray_features_standard_vp.rds',
  microarray.modules = snakemake@input[[14]],#'Data/Microarray/analysis_output/stability/microarray_modules_standard_vp.rds',
  tbnks = snakemake@input[[15]]#'Data/TBNK/analysis_output/stability/tbnks_features_standard_vp.rds'
)

## Sample level esets, subset to just stable features
STABLE.ESET.SAMPLE.LEVEL.OUT.PATHS = list(
  somalogic.features = snakemake@output[[1]],#'Data/Somalogic/analysis_output/stability/stable_somalogic_sample_level_features.rds',
  somalogic.modules = snakemake@output[[2]],#'Data/Somalogic/analysis_output/stability/stable_somalogic_sample_level_modules.rds',
  microarray.features = snakemake@output[[3]],#'Data/Microarray/analysis_output/stability/stable_microarray_sample_level_features.rds',
  microarray.modules = snakemake@output[[4]],#'Data/Microarray/analysis_output/stability/stable_microarray_sample_level_modules.rds',
  tbnks = snakemake@output[[5]]#'Data/TBNK/analysis_output/stability/stable_tbnk_sample_level_features.rds'
)

## Subject level esets, subset to just stable features
STABLE.ESET.SUBJECT.LEVEL.OUT.PATHS = list(
  somalogic.features = snakemake@output[[6]],#'Data/Somalogic/analysis_output/stability/stable_somalogic_subject_level_features.rds',
  somalogic.modules = snakemake@output[[7]],#'Data/Somalogic/analysis_output/stability/stable_somalogic_subject_level_modules.rds',
  microarray.features = snakemake@output[[8]],#'Data/Microarray/analysis_output/stability/stable_microarray_subject_level_features.rds',
  microarray.modules = snakemake@output[[9]],#'Data/Microarray/analysis_output/stability/stable_microarray_subject_level_modules.rds',
  tbnks = snakemake@output[[10]]#'Data/TBNK/analysis_output/stability/stable_tbnk_subject_level_features.rds'
)

# Load the data
sample.esets = lapply(ESET.SAMPLE.LEVEL.IN.PATHS, readRDS)
subject.esets = lapply(ESET.SUBJECT.LEVEL.IN.PATHS, readRDS)
vps = lapply(VP.IN.PATHS, readRDS)

# Make a function to subset an eset to just the features for which
# the patient covariate explains at least half the total variance
select.stable = function(eset, vp) {
  stable.eset = eset[rownames(vp)[vp$Patient >= .5], ]
  return(stable.eset)
}

# Initate a function to save the esets
save.eset = function(eset, path) {
  saveRDS(eset, path)
  return('eset saved')
}

# Subset esets
stable.sample.esets = mapply(select.stable, sample.esets, vps, SIMPLIFY = FALSE)
stable.subject.esets = mapply(select.stable, subject.esets, vps, SIMPLIFY = FALSE)

# Save esets
out = mapply(save.eset, stable.sample.esets, STABLE.ESET.SAMPLE.LEVEL.OUT.PATHS)
out = mapply(save.eset, stable.subject.esets, STABLE.ESET.SUBJECT.LEVEL.OUT.PATHS)