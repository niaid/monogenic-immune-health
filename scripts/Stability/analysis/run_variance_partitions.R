## Run the patient-residual variance partitions to estimate feature stability

# Set Seed
set.seed(131)

# Load libraries
library(variancePartition)
library(Biobase)
source('scripts/util/VariancePartition/variancePartition.R')

# Set input paths
## Sample-level eset input paths
ESET.IN.PATHS = list(
  somalogic.features = snakemake@input[[1]],#'Data/Somalogic/data_analysis_ready/analysis_ready_sample_level_training_somalogic.rds',
  somalogic.modules = snakemake@input[[2]],#'Data/Somalogic/analysis_output/wgcna_results/scores_sample_level.rds',
  microarray.features = snakemake@input[[3]],#'Data/Microarray/data_analysis_ready/eset_batch_training_sample.rds',
  microarray.modules = snakemake@input[[4]],#'Data/Microarray/analysis_output/WGCNA/array_sample_scores.rds',
  tbnks = snakemake@input[[5]]#'Data/TBNK/data_analysis_ready/tbnk_training_sample_level_eset.rds'
)

# Set output paths
## Variance partitions into patient-explained variance and residual variance
VP.OUT.PATHS = list(
  somalogic.features = snakemake@output[[1]],#'Data/Somalogic/analysis_output/stability/somalogic_features_standard_vp.rds',
  somalogic.modules = snakemake@output[[2]],#'Data/Somalogic/analysis_output/stability/somalogic_modules_standard_vp.rds',
  microarray.features = snakemake@output[[3]],#'Data/Microarray/analysis_output/stability/microarray_features_standard_vp.rds',
  microarray.modules = snakemake@output[[4]],#'Data/Microarray/analysis_output/stability/microarray_modules_standard_vp.rds',
  tbnks = snakemake@output[[5]]#'Data/TBNK/analysis_output/stability/tbnks_features_standard_vp.rds'
)

# Load data
esets = lapply(ESET.IN.PATHS, readRDS)

# Run variance partitions
print(VP.OUT.PATHS)
vps = lapply(esets, patient_id_variance_partition)

# Save results
mapply(function(vp, out.path) {
  saveRDS(vp, out.path)
}, vps, VP.OUT.PATHS)