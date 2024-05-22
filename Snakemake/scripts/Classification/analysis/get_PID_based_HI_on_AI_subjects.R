# Derive the HI scores for PID subjects based on the classifier that was only trained on AI and healthy

# Set globals
## Random Forest models trained on PID and Healthy Patients
PID.MODELS.PATH = snakemake@input[[1]] #'Classification/results/healthy_rf_models_PID.RDS'
## Random Forest design matrices for all subjects
FULL.DESGIN.MATRICES.PATH = snakemake@input[[2]] #'Classification/design_matrices/healthy_random_forest_design_matrices_all.RDS'
## Random Forest sample meta data for all subjects
FULL.META.DATA.PATH = snakemake@input[[3]] #'Classification/design_matrices/healthy_random_forest_sample_meta_data_all.RDS'

## The predictive index of each subject
SCORES.OUT.PATH = snakemake@output[[1]] #'Classification/predictions/healthy_rf_AI_predictions_using_PID_index.RDS'

# Load libraries
library(randomForest)
library(Biobase)

# Set seed
set.seed(102409)

# Source utility functions
source('scripts/util/Groups/groups.R')

# Load data
models = readRDS(PID.MODELS.PATH)

full.design.matrices = readRDS(FULL.DESGIN.MATRICES.PATH)
full.meta.data = readRDS(FULL.META.DATA.PATH)

# For each model
scores = mapply(function(model, full.design.matrix) {
  # Make the design matrix for just the AI patients
  X = full.design.matrix[full.meta.data$condition %in% util.get_ai(), ]
  
  # Get the classifier's predictions
  predict(model, X, type = 'prob')[,'1']
}, models, full.design.matrices)

# Save results
saveRDS(as.data.frame(scores), SCORES.OUT.PATH)