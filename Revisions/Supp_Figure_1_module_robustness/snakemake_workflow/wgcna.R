## Create microarray WGCNA modules and save the variances explained by
## PC1 of each module

## Load Libraries
library(WGCNA)
library(Biobase)

boot_iter <- as.numeric(snakemake@wildcards$boot)
set.seed(boot_iter)

## Source wgcna function
source('../../../2021_07_08/scripts/util/WGCNA/runWGCNA.r')
source('../../../2021_07_08/scripts/util/Processing/averageRepeatSamples.R')
source('../../../2021_07_08/scripts/util/WGCNA/get_eigengene_scores.R')
source('../../../2021_07_08/scripts/util/Processing/removeOutlierPatients.R')

# Set GlobalVariables
## Clean sample-level somalogic data
SAMPLES.IN.PATH = snakemake@input[[1]]#'Data/Microarray/data_analysis_ready/eset_batch_training_sample.rds'

## The somalogic WGCNA feature to module map
MODULES.OUT.PATH = snakemake@output[[1]]#'Data/Microarray/analysis_output/WGCNA/modules.rds'
## The sample-level somalogic module scores
SCORES.SAMPLE.LEVEL.OUT.PATH = snakemake@output[[2]]#'Data/Microarray/analysis_output/WGCNA/array_sample_scores.rds'
## The subject-level somalogic module scores
SCORES.SUBJECT.LEVEL.OUT.PATH = snakemake@output[[3]]#'Data/Microarray/analysis_output/WGCNA/array_subject_scores.rds'
## The variances explained by PC1 of each module

## Load data
microarray.samples = readRDS(SAMPLES.IN.PATH)

## Prevent WGCNA from operating with parallel
disableWGCNAThreads()

## Remove outlier samples (and plot results)
microarray.samples.filtered = removeOutlierPatients(microarray.samples, cutHeight = 250)

## Calculate the subject level data without outliers
microarray.subjects = averageRepeatSamples(microarray.samples.filtered)


keep_subj <- sample.int(ncol(microarray.subjects), round(.8 * ncol(microarray.subjects)))
microarray.subjects = microarray.subjects[, keep_subj]

## Run wgcna function
if(snakemake@wildcards$datatype == "array"){
  modules = runWGCNA(microarray.subjects, OUTDIR, method = 'hybrid', pamStage = TRUE, 
                     pamRespectsDendro = FALSE, beta = 12, minModuleSize = 30, deepSplit = 2,
                     intermediate.results.path = NULL, 
                     diagnostic.plots.path = NULL)

}else if(snakemake@wildcards$datatype == "soma"){
  modules = runWGCNA(microarray.subjects, OUTDIR, method = 'tree', pamStage = FALSE, pamRespectsDendro = FALSE,
                   beta = 12, deepSplit = 2, minModuleSize = 30, 
                   intermediate.results.path = NULL, 
                   diagnostic.plots.path = NULL)
}

## Get the scores associated with each sample for each module
scores.sample.level = get_eigengene_scores(microarray.subjects, microarray.samples, modules)

## Get the variances explained by PC1 of each module
variances = get_eigengene_variance_explained(microarray.subjects, modules)

## Average over repeat samples
scores.subject.level = averageRepeatSamples(scores.sample.level)

## Save modules and scores
saveRDS(modules, file = MODULES.OUT.PATH)
saveRDS(scores.sample.level, file = SCORES.SAMPLE.LEVEL.OUT.PATH)
saveRDS(scores.subject.level, file = SCORES.SUBJECT.LEVEL.OUT.PATH)
#saveRDS(variances, file = VARIANCES.OUT.PATH)
