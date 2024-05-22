# Computes the somalogic module scores for the testing set

# Source utility functions
source('scripts/util/Processing/averageRepeatSamples.R')
source('scripts/util/WGCNA/get_eigengene_scores.R')

# Set paths
## Microarray modules
MODULES.IN.PATH = snakemake@input[[1]]#'Data/Microarray/analysis_output/WGCNA/modules.rds'
## Microarray subject-level training eset
TRAINING.SET.MICROARRAY.ESET.IN.PATH = snakemake@input[[2]]#'Data/Microarray/data_analysis_ready/eset_batch_training_subject.rds'
## Microarray sample-level testing eset
TESTING.SET.MICROARRAY.ESET.IN.PATH = snakemake@input[[3]]#'Data/Microarray/data_analysis_ready/eset_batch_validation_sample.rds'

TESTING.SET.SAMPLE.LEVEL.SCORES.ESET.OUT.PATH = snakemake@output[[1]]#'Data/Microarray/analysis_output/WGCNA/scores_sample_level_testing.rds'
TESTING.SET.SUBJECT.LEVEL.SCORES.ESET.OUT.PATH = snakemake@output[[2]]#'Data/Microarray/analysis_output/WGCNA/scores_subject_level_testing.rds'

# Load data
modules = readRDS(MODULES.IN.PATH)
training.set.microarray.eset = readRDS(TRAINING.SET.MICROARRAY.ESET.IN.PATH)
testing.set.microarray.eset = readRDS(TESTING.SET.MICROARRAY.ESET.IN.PATH)

# Get the sample level module scores for the testing eset
testing.set.sample.level.scores.eset = get_eigengene_scores(training.set.microarray.eset, testing.set.microarray.eset, modules)

# Average over samples within a subject to get the subject level module scores for the testing eset
testing.set.subject.level.scores.eset = averageRepeatSamples(testing.set.sample.level.scores.eset)

# Save results
saveRDS(testing.set.sample.level.scores.eset, TESTING.SET.SAMPLE.LEVEL.SCORES.ESET.OUT.PATH)
saveRDS(testing.set.subject.level.scores.eset, TESTING.SET.SUBJECT.LEVEL.SCORES.ESET.OUT.PATH)
