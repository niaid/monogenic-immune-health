# Computes the somalogic module scores for the testing set

# Source utility functions
source('scripts/util/Processing/averageRepeatSamples.R')
source('scripts/util/WGCNA/get_eigengene_scores.R')

# Set paths
## Somalogic modules
MODULES.IN.PATH = snakemake@input[[1]]#'Data/Somalogic/analysis_output/wgcna_results/modules.rds'
## Somalogic subject-level training eset
TRAINING.SET.SOMALOGIC.ESET.IN.PATH = snakemake@input[[2]]#'Data/Somalogic/data_analysis_ready/analysis_ready_subject_level_training_somalogic.rds'
## Somalogic sample-level testing eset
TESTING.SET.SOMALOGIC.ESET.IN.PATH = snakemake@input[[3]]#'Data/Somalogic/data_analysis_ready/analysis_ready_sample_level_testing_somalogic.rds'

## Sample level somalogic module scores for testing set
TESTING.SET.SAMPLE.LEVEL.SCORES.ESET.OUT.PATH = snakemake@output[[1]]#'Data/Somalogic/analysis_output/wgcna_results/scores_sample_level_testing.rds'
## Subject level somalogic module scores for testing set
TESTING.SET.SUBJECT.LEVEL.SCORES.ESET.OUT.PATH = snakemake@output[[2]]#'Data/Somalogic/analysis_output/wgcna_results/scores_subject_level_testing.rds'

# Load data
modules = readRDS(MODULES.IN.PATH)
training.set.somalogic.eset = readRDS(TRAINING.SET.SOMALOGIC.ESET.IN.PATH)
testing.set.somalogic.eset = readRDS(TESTING.SET.SOMALOGIC.ESET.IN.PATH)

# Get the sample level module scores for the testing eset
testing.set.sample.level.scores.eset = get_eigengene_scores(training.set.somalogic.eset, testing.set.somalogic.eset, modules)

# Average over samples within a subject to get the subject level module scores for the testing eset
testing.set.subject.level.scores.eset = averageRepeatSamples(testing.set.sample.level.scores.eset)

# Save results
saveRDS(testing.set.sample.level.scores.eset, TESTING.SET.SAMPLE.LEVEL.SCORES.ESET.OUT.PATH)
saveRDS(testing.set.subject.level.scores.eset, TESTING.SET.SUBJECT.LEVEL.SCORES.ESET.OUT.PATH)