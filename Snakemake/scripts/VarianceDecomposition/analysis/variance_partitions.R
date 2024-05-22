## Run variance partition models with patients, conditions, and medications
## as covariates for each data type

# Load libraries
library(ggplot2)
library(limma)
library(doParallel)
library(variancePartition)
library(Biobase)

# Source utilities
source('scripts/util/VariancePartition/variancePartition.R')

# Set Globals
## Path to the medication inforamtion on patients
MEDICATIONS.IN.PATH = snakemake@input[[1]]#'Medications/medications.types.rds'

## Path to the sample-level expression sets for each dataset
ESETS.IN.PATHS = list(
  somalogic.features = snakemake@input[[2]],#'Data/Somalogic/data_analysis_ready/analysis_ready_sample_level_training_somalogic.rds',
  somalogic.modules = snakemake@input[[3]],#'Data/Somalogic/analysis_output/wgcna_results/scores_sample_level.rds',
  microarray.features = snakemake@input[[4]],#'Data/Microarray/data_analysis_ready/eset_batch_training_sample.rds',
  microarray.modules = snakemake@input[[5]],#'Data/Microarray/analysis_output/WGCNA/array_sample_scores.rds',
  tbnks = snakemake@input[[6]]#'Data/TBNK/data_analysis_ready/tbnk_training_sample_level_eset.rds'
)

## Paths to save the variance partition results
VP.OUT.PATHS = list(
  somalogic.features = snakemake@output[[1]],#'Data/Somalogic/analysis_output/variance_decomposition/somalogic_features_vp.RDS',
  somalogic.modules = snakemake@output[[2]],#'Data/Somalogic/analysis_output/variance_decomposition/somalogic_modules_vp.RDS',
  microarray.features = snakemake@output[[3]],#'Data/Microarray/analysis_output/variance_decomposition/microarray_features_vp.RDS',
  microarray.modules = snakemake@output[[4]],#'Data/Microarray/analysis_output/variance_decomposition/microarray_modules_vp.RDS',
  tbnks = snakemake@output[[5]]#'Data/TBNK/analysis_output/variance_decomposition/tbnk_features_vp.RDS'
)

# Load data
medications = readRDS(MEDICATIONS.IN.PATH)
esets = lapply(ESETS.IN.PATHS, readRDS)

# Remove unwanted medications and patient info columns from the medications matrix
medications.to.remove = c('Procedure', 'Surgery', 'Transfusion', 'Transplantation', 'Other')
columns.to.remove = c('patient_id','visit_id', medications.to.remove)
medications = medications[, setdiff(colnames(medications), columns.to.remove)]

# Get a matrix of medications that corresponds to the visits in each data set
medication.matrices = lapply(esets, function(eset) {
  medications[colnames(eset), ]
})

# Run variance partitions
vps = mapply(function(eset, medication.matrix) {
  condition_medication_variance_partition(eset, medication.matrix)
}, esets, medication.matrices, SIMPLIFY = FALSE)

# Save results
mapply(function(vp, out.path) {saveRDS(vp, out.path)}, vps, VP.OUT.PATHS)