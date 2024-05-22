## Create microarray WGCNA modules and save the variances explained by
## PC1 of each module

## Load Libraries
library(WGCNA)
library(Biobase)

## Set seed
set.seed(10798)

## Source wgcna function
source('scripts/util/WGCNA/runWGCNA.r')
source('scripts/util/Processing/averageRepeatSamples.R')
source('scripts/util/WGCNA/get_eigengene_scores.R')
source('scripts/util/Processing/removeOutlierPatients.R')

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
VARIANCES.OUT.PATH = snakemake@output[[4]]#'Data/Microarray/analysis_output/WGCNA/array_subject_varexp.rds'
## WGCNA intermediate objects
INTERMEDIATES.OUT.PATH = snakemake@output[[5]]#'Data/Microarray/analysis_output/WGCNA/WGCNA_microarray_intermediates.rds'
## Outlier removal plots
OUTLIER.REMOVAL.PLOTS.OUT.PATH = snakemake@output[[6]]#'Paper_1_Figures/Supplemental_Figure_1/microarray_outlier_removal_for_wgcna.pdf'
## Diagnostic plots from WGCNA module creation
WGCNA.PLOTS.OUT.PATH = snakemake@output[[7]]#'Paper_1_Figures/Supplemental_Figure_1/microarray_wgcna_module_creation.pdf'

## Load data
microarray.samples = readRDS(SAMPLES.IN.PATH)

## Prevent WGCNA from operating with parallel
disableWGCNAThreads()

## Remove outlier samples (and plot results)
pdf(OUTLIER.REMOVAL.PLOTS.OUT.PATH)
microarray.samples.filtered = removeOutlierPatients(microarray.samples, cutHeight = 250)
dev.off()

## Calculate the subject level data without outliers
microarray.subjects = averageRepeatSamples(microarray.samples.filtered)

## Run wgcna function
modules = runWGCNA(microarray.subjects, OUTDIR, method = 'hybrid', pamStage = TRUE, 
                   pamRespectsDendro = FALSE, beta = 12, minModuleSize = 30, deepSplit = 2,
                   intermediate.results.path = INTERMEDIATES.OUT.PATH, 
                   diagnostic.plots.path = WGCNA.PLOTS.OUT.PATH)

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
saveRDS(variances, file = VARIANCES.OUT.PATH)
