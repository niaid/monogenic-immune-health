## Predicts the HI scores for the testing set based on previously trained RF models.
## The foreground condition group (elsewhere just called the "condition group") 
## and the background condition group are specified via group keys in the file name in the snakemake pipeline.
## Group keys are mapped to vectors of conditions in the script create_classifier_groups.R. 
## Only output from the corresponding RF analysis this condition group and background condition group 
## are included in the output from this script.

# Set globals
## Random Forest models trained on training subjects
TRAINING.MODELS.PATH = snakemake@input[[1]] #'Classification/healthy_random_forest_design_matrices_all.RDS'
## Random Forest design matrices for testing subjects
TESTING.DESGIN.MATRICES.PATH = snakemake@input[[2]] #'Classification/healthy_random_forest_testing_design_matrices_all.RDS'

## The predictive index of each subject
SCORES.OUT.PATH = snakemake@output[[1]] #'Classification/predictions/healthy_rf_testing_predictions_all.RDS'

# Load libraries
library(randomForest)
library(Biobase)

# Set seed
set.seed(102409)

# Source utility functions
source('scripts/util/Groups/groups.R')

# Load data
models = readRDS(TRAINING.MODELS.PATH)
testing.design.matrices = readRDS(TESTING.DESGIN.MATRICES.PATH)

# For each model
scores = mapply(function(model, testing.design.matrix) {
  # Get the classifier's predictions
  predict(model, testing.design.matrix, type = 'prob')[,'1']
}, models, testing.design.matrices)

# Save results
saveRDS(as.data.frame(scores), SCORES.OUT.PATH)