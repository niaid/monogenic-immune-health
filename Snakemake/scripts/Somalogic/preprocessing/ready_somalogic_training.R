## Remove outliers / unwanted samples and
## unwanted somamers from the somalogic training data
## and also create another data set with repeat samples
## averaged within each subject

# Load Libraries
library(WGCNA)
library(Biobase)

# Source utility functions
source('scripts/util/Processing/averageTechnicalReplicates.R')
source('scripts/util/Processing/averageRepeatSamples.R')

# Set Global Variables
## The somalogic training eset (prior to cleaning)
SAMPLES.IN.PATH = snakemake@input[[1]]#'Data/Somalogic/processed/v1/training_somalogic.rds'

## The cleaned up sample-level somalogic training eset
SAMPLE.LEVEL.OUT.PATH = snakemake@output[[1]]#'Data/Somalogic/data_analysis_ready/analysis_ready_sample_level_training_somalogic.rds'
## The cleaned up sample-level somalogic training eset
SUBJECT.LEVEL.OUT.PATH = snakemake@output[[2]]#'Data/Somalogic/data_analysis_ready/analysis_ready_subject_level_training_somalogic.rds'

# Load data
somalogic.samples = readRDS(SAMPLES.IN.PATH)

# Remove unwanted visits for the following reasons:
## V313: hemolysis level 4
## V282: hemolysis level 4
## V210: very odd cloudy/milky sample
samples.to.remove = c('V313', 'V282', 'V210')
somalogic.samples = somalogic.samples[, ! somalogic.samples$visit_id %in% samples.to.remove]

# Remove unwanted somamers for the following reasons
## EGFRvIII: Removed from Somalogic panel, found to be cross reactive to an unknown source
features.to.remove = c('EGFRvIII')
somalogic.samples = somalogic.samples[! rownames(somalogic.samples) %in% features.to.remove,]

# Average technical replicates
somalogic.samples = averageTechnicalReplicates(somalogic.samples,
                                               visit.id.col = 'visit_id',
                                               meta.cols = c('patient_id',
                                                             'gender',
                                                             'patient_age_at_time_of_blood_draw',
                                                             'race',
                                                             'condition',
                                                             'ethnicity',
                                                             'plate_id',
                                                             'assay_desc',
                                                             'visit_type'))

# Rename columns to the visit id
colnames(somalogic.samples) = somalogic.samples$visit_id

# Collapse samples into subject through averaging
somalogic.subjects = averageRepeatSamples(somalogic.samples)

# Save results
saveRDS(somalogic.samples, file = SAMPLE.LEVEL.OUT.PATH)
saveRDS(somalogic.subjects, file = SUBJECT.LEVEL.OUT.PATH)