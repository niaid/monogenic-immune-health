## This script uses regularized logistic regression to create a signature of the most important
## genes for classification of healthy versus non-healthy subjects. We also get a CV AUC estimate to make
## sure the classifier is reasonable on the training data. This signature is used a positive control to ensure we aren’t losing 
## too much information when transferring from the various data types to transcriptional data.

# Set seed
set.seed(140)

# Load libraries and source utility functions
library(Biobase)
library(glmnet)
source('scripts/util/Classification/get_aucs.R')

# Set paths
## Design matrices with various data types
MATRICES.IN.PATH = snakemake@input[[1]]#'Classification/design_matrices/healthy_all_design_matrices_all.RDS'
## Corresponding metdata
META.IN.PATH = snakemake@input[[2]]#'Classification/meta_data/healthy_random_forest_sample_meta_data_all.RDS'

## Microarray classifier gene signature
SIGNATURE.OUT.PATH = snakemake@output[[1]]#'Classification/transcriptional_surrogates/microarray_classifier_signatures.RDS'

# Load design matrices
Xs = readRDS(MATRICES.IN.PATH)
meta = readRDS(META.IN.PATH)

# Extract microarray design matrix 
X = Xs$microarray.features

# Get response vector
y = as.numeric(meta$condition == 'Healthy')

# Scale design matrix
X = scale(X)

# Divide indices into n randomly sampled cross validation groups
xs = 1:length(y)
n = 10
f = factor(xs %% n)
groups = split(sample(xs), f)

# Her we get the cross validation predictions
# For each group of samples
predictions = lapply(groups, function(group) {
  # Subset the design matrix to samples not in that group for training
  X.train = X[-group, , drop = FALSE]
  # Subset the design matrix to just samples in that group for testing
  X.test = X[group, , drop = FALSE]
  # Get the responses vector of the training group
  y.train = y[-group]
  # Train a L2 penalty logistic regression model using the training data
  model = cv.glmnet(X.train, y.train, family = 'binomial', alpha = 0)
  # Get the probabilities that each sample from the test group is of the positive class (healthy)
  predictions = predict(model, X.test)[,'1']
})

# Put predictions into a vector and get corresponding conditions
predictions = unname(unlist(predictions))
conditions = meta$condition[unlist(groups)]

# Get an estimate for the auc
roc = get_roc(predictions, conditions, 'Healthy')
auc = get_auc(roc)
print('Classifier AUC:')
print(auc)

# Get most important m features for the prediction
m = 500
model = cv.glmnet(X, y, family = 'binomial', alpha = 0)

# Here we extract the model coefficients (note that the coefficients with
# highest absolute values should be the most important, as the data was scaled)
coefs = coef(model)
# Convert the coefficents to a vector
coefs = coefs[,1]
# Sort the coefs by absolute value
coefs = coefs[order(abs(coefs), decreasing = TRUE)]
# Remove the intercept term
coefs = coefs[names(coefs) != '(Intercept)']
# Remove the 'microarray.features' prefixes from the genes
names(coefs) = gsub('microarray\\.features\\.', '', names(coefs))
# Extract the top m features from the model
features = names(coefs)[1:m]

# Make a signature based on the most important features
signature = list(positive = features[coefs[features] > 0], negative = features[coefs[features] < 0])

# Save the signature
saveRDS(signature, SIGNATURE.OUT.PATH)
