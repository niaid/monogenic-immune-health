## This script takes an input and generates permutations of the
## design matrices used to obtain the GVIs for the random forest classifiers.
## p-values are then assigned to each feature for each classifier based on the
## permutation testing.
## The foreground condition group (elsewhere just called the "condition group") 
## and the background condition group are specified via group keys in the file name in the snakemake pipeline.
## Group keys are mapped to vectors of conditions in the script create_classifier_groups.R. 
## Only output from the corresponding RF analysis this condition group and background condition group 
## are included in the output from this script.

# Set seed
set.seed(99)

# Load packages and source utilities
library(randomForest)
source('scripts/util/Classification/randomForestClassifier.R')

# Set globals
DESIGN.MATRICES.IN.PATH = snakemake@input[[1]]#'Classification/healthy_random_forest_design_matrices_all.RDS'
META.DATA.IN.PATH = snakemake@input[[2]]#'Classification/healthy_random_forest_sample_meta_data_all.RDS'
GVIS.IN.PATH = snakemake@input[[3]]#'Classification/results/healthy_random_forest_rf_gvis_all.RDS'
CONDITION.GROUPS.IN.PATH = snakemake@input[[4]]#'Classification/condition_groups.RDS'
BACKGROUND.GROUPS.IN.PATH = snakemake@input[[5]]#'Classification/background_groups.RDS'

SNAKEMAKE.OUT.PATH = snakemake@output[[1]]

# We get the group from the output file name
out.file = basename(SNAKEMAKE.OUT.PATH)
out.file = gsub('.RDS$', '', out.file, ignore.case = T)
fields = strsplit(out.file,'_')[[1]]
condition.id = fields[1] # The condition group represents the 'positive' condition for the classifier (e.g. 'healthy', 'cgd', 'xcgd','47cgd','stat1.gof')
background.id = fields[length(fields)] # The background group represents the background pool of all other conditions to consider (e.g. 'PID','AI','all')
i = as.numeric(fields[3])

# We load the condition and background groups
condition.groups = readRDS(CONDITION.GROUPS.IN.PATH)
background.groups = readRDS(BACKGROUND.GROUPS.IN.PATH)

# Get the seed corresponding to each condition so we aren't using the same seeds for each condition
n = 100000
condition.seeds = lapply(condition.groups, function(condition.group){sample.int(n, size = 1)})
background.seeds = lapply(background.groups, function(background.group){sample.int(n, size = 1)})

condition.seed = condition.seeds[[condition.id]]
background.seed = background.seeds[[background.id]]

# Get the conditons to investigate corresponding to the 'condition' field.
condition = condition.groups[[condition.id]]

# Set seed based on that permutation number
set.seed(background.seed + condition.seed + i)

# We load the design matrices, associated meta data, and gvis
# from the random forest classifiers
Xs = readRDS(DESIGN.MATRICES.IN.PATH)
meta = readRDS(META.DATA.IN.PATH)
gvis = readRDS(GVIS.IN.PATH)

# For each design matrix and its associated RF GVIs
perm.pvals = mapply(function(X, gvi) {
  # Repeat for the following for the number of features
  # in the design matrix:
  perm.gvis = sapply(1:ncol(X), function(j) {
    # Permute the response vector
    y.perm = sample(meta$condition)
    # Train a RF model using the permuted
    # response vector, and return the associated gvis
    # for each feature
    get.gvis(X, y.perm, pos = condition)
  })
  
  # For each feature, get the permutation pvalue, which is the 
  # percent of times when the true gvi for that feature was 
  # lower than the permutation gvis for 
  # that feature in this
  # iteration of the permutation tests
  rowMeans(gvi <= perm.gvis)

}, Xs, gvis, SIMPLIFY = FALSE)

# Save results in the output file for
# this permutation number, as specified
# in the snakefile
saveRDS(perm.pvals, SNAKEMAKE.OUT.PATH)

# Use the permutation number to print 
# confirmation of finishing this permutation
print(i)