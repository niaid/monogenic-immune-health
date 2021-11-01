## Creates somalogic WGCNA modules. Calculates the module PC1 score of each sample
## (with the PC loadings themselves derived at the subject level as to not bias toward patients
## with multiple samples). Also averages over these sample-level scores to create a subject-level
## scores eset. Also saves the variances explained by PC1 of each module (derived at the subject-level).

# Load Libraries
library(WGCNA)
library(Biobase)

# Set seed
set.seed(130)

# Source wgcna utility functions
source('scripts/util/WGCNA/runWGCNA.r')
source('scripts/util/Processing/averageRepeatSamples.R')
source('scripts/util/WGCNA/get_eigengene_scores.R')
source('scripts/util/Processing/removeOutlierPatients.R')

# Set GlobalVariables
## Clean sample-level somalogic data
SAMPLES.IN.PATH = snakemake@input[[1]]#'Data/Somalogic/data_analysis_ready/analysis_ready_sample_level_training_somalogic.rds'

## The somalogic WGCNA feature to module map
MODULES.OUT.PATH = snakemake@output[[1]]#'Data/Somalogic/analysis_output/wgcna_results/modules.rds'
## The sample-level somalogic module scores
SCORES.SAMPLE.LEVEL.OUT.PATH = snakemake@output[[2]]#'Data/Somalogic/analysis_output/wgcna_results/scores_sample_level.rds'
## The subject-level somalogic module scores
SCORES.SUBJECT.LEVEL.OUT.PATH = snakemake@output[[3]]#'Data/Somalogic/analysis_output/wgcna_results/scores_subject_level.rds'
## The variances explained by PC1 of each module
VARIANCES.OUT.PATH = snakemake@output[[4]]#'Data/Somalogic/analysis_output/wgcna_results/variances.rds'
## WGCNA intermediate objects
INTERMEDIATES.OUT.PATH = snakemake@output[[5]]#'Data/Somalogic/analysis_output/wgcna_results/WGCNA_somalogic_intermediates.rds'
## Outlier removal plots
OUTLIER.REMOVAL.PLOTS.OUT.PATH = snakemake@output[[6]]#'Paper_1_Figures/Supplemental_Figure_1/somalogic_outlier_removal_for_wgcna.pdf'
## Diagnostic plots from WGCNA module creation
WGCNA.PLOTS.OUT.PATH = snakemake@output[[7]]#'Paper_1_Figures/Supplemental_Figure_1/somalogic_wgcna_module_creation.pdf'

# Load data
somalogic.samples = readRDS(SAMPLES.IN.PATH)

# Prevent WGCNA from operating with parallel
disableWGCNAThreads()

# Remove outlier samples (and plot results)
pdf(OUTLIER.REMOVAL.PLOTS.OUT.PATH)
somalogic.samples.filtered = removeOutlierPatients(somalogic.samples, cutHeight = 75)
dev.off()

# Calculate the subject level data without outliers
somalogic.subjects = averageRepeatSamples(somalogic.samples.filtered)

# Run wgcna function
## Here we use pamStage = FALSE and method = 'tree' even though these options are set to TRUE
## and 'hybrid' in the tutorial respectively,
## because using the tutorial's options result in lower median module variances explained and a smaller
## number of modules
modules = runWGCNA(somalogic.subjects, OUTDIR, method = 'tree', pamStage = FALSE, pamRespectsDendro = FALSE,
                   beta = 12, deepSplit = 2, minModuleSize = 30, 
                   intermediate.results.path = INTERMEDIATES.OUT.PATH, 
                   diagnostic.plots.path = WGCNA.PLOTS.OUT.PATH)

# Get the scores associated with each sample for each module
scores.sample.level = get_eigengene_scores(somalogic.subjects, somalogic.samples, modules)

# Get the variances explained by PC1 of each module
variances = get_eigengene_variance_explained(somalogic.subjects, modules)

# Average over repeat samples
scores.subject.level = averageRepeatSamples(scores.sample.level)

# Save modules and scores
saveRDS(modules, file = MODULES.OUT.PATH)
saveRDS(scores.sample.level, file = SCORES.SAMPLE.LEVEL.OUT.PATH)
saveRDS(scores.subject.level, file = SCORES.SUBJECT.LEVEL.OUT.PATH)
saveRDS(variances, file = VARIANCES.OUT.PATH)

