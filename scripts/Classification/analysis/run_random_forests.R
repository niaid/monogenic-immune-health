## This gets the LOO CV predictions of a subject being in one of the specified foreground groups based on
## a random forest classifier, for each of the selected design matrices.
## Additionally, for each design matrix, this script runs a random forest 
## classifier using all samples, and gets the resulting GVIs
## associated with each feature.
## It also saves the trained RF model from which the GVIs are derived.
## The foreground condition group (elsewhere just called the "condition group") 
## and the background condition group are specified via group keys in the file name in the snakemake pipeline.
## Group keys are mapped to vectors of conditions in the script create_classifier_groups.R. 
## Only conditions in the condition group and background condition group are included in the output from this script.

# Set seed
set.seed(2020)

# Load packages and source utilities
library(randomForest)
source('scripts/util/Classification/randomForestClassifier.R')

# Set globals
DESIGN.MATRICES.IN.PATH = snakemake@input[[1]]#'Classification/healthy_random_forest_design_matrices_all.RDS'
META.DATA.IN.PATH = snakemake@input[[2]]#'Classification/healthy_random_forest_sample_meta_data_all.RDS'
CONDITION.GROUPS.IN.PATH = snakemake@input[[3]]#Classification/condition_groups.RDS
BACKGROUND.GROUPS.IN.PATH = snakemake@input[[4]]#Classification/background_groups.RDS
PREDICTIONS.OUT.PATH = snakemake@output[[1]]#'Classification/results/healthy_rf_results_all.RDS'
GVIS.OUT.PATH = snakemake@output[[2]]#'Classification/results/healthy_rf_gvis_all.RDS'
MODELS.OUT.PATH = snakemake@output[[3]]#'Classification/results/healthy_rf_models_all.RDS'

# We get the group from the output file name
out.file = basename(DESIGN.MATRICES.IN.PATH)
out.file = gsub('.RDS$', '', out.file, ignore.case = T)
fields = strsplit(out.file,'_')[[1]]
condition.id = fields[1] # The condition group represents the 'positive' condition for the classifier (e.g. 'healthy', 'cgd', 'xcgd','47cgd','stat1.gof')
background.id = fields[length(fields)] # The background group represents the background pool of all other conditions to consider (e.g. 'PID','AI','all')

# We load the condition and background groups
condition.groups = readRDS(CONDITION.GROUPS.IN.PATH)
background.groups = readRDS(BACKGROUND.GROUPS.IN.PATH)

# Get the seed corresponding to each condition so we aren't using the same seeds for each condition
n = 100000
condition.seeds = lapply(condition.groups, function(condition.group){sample.int(n, size = 1)})
background.seeds = lapply(background.groups, function(background.group){sample.int(n, size = 1)})

condition.seed = condition.seeds[[condition.id]]
background.seed = background.seeds[[background.id]]

# Get the conditions to serve as the 'positive' class
condition.group = condition.groups[[condition.id]]

# Set the first seed
set.seed(sample.int(n, size = 1) + condition.seed + background.seed)

# Load data
Xs = readRDS(DESIGN.MATRICES.IN.PATH)
meta = readRDS(META.DATA.IN.PATH)

# Run cross validations on each design matrix
predictions = sapply(Xs, function(X) {
  cross.validation(X, meta$condition, pos = condition.group)
})

# Name the columns and rows of the predictions matrix
colnames(predictions) = names(Xs)
rownames(predictions) = rownames(meta)

# Convert the results matrix to a data frame
predictions = as.data.frame(predictions)

# Reset seed
set.seed(sample.int(n, size = 1) + condition.seed + background.seed)

# Train the random forests with all models
models = lapply(Xs, function(X) {
  get.rf.model(X, meta$condition, pos = condition.group)
})

# Get the gvis from random forests with all samples
gvis = lapply(models, function(model) {
  model$importance[,'MeanDecreaseGini']
})

# Save the results
saveRDS(predictions, PREDICTIONS.OUT.PATH)
saveRDS(gvis, GVIS.OUT.PATH)
saveRDS(models, MODELS.OUT.PATH)
