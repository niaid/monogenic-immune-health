## Creates microarray surrogate signatures for a large number of non-microarray features,
## namely the HI driving features, the HI itself, and the jPC1. We also include pseudo-signatures
## (sets of genes in the same format as the signatures, but that actually came from microarray data,
## and thus don't require conversion) such as the important microarray modules, and the microarray-based
## logistic regression classifier (the "positive control")

# Set globals
## Design matrices
MATRICES.IN.PATH = snakemake@input[[1]]#'"Classification/design_matrices/healthy_all_design_matrices_all.RDS"'
## Corresponding meta data
META.IN.PATH = snakemake@input[[2]]#'Classification/meta_data/healthy_random_forest_sample_meta_data_all.RDS'
## GVI permutation testing pvalues
PVALS.IN.PATH = snakemake@input[[3]]#'Classification/results/healthy_rf_pvals_all.RDS'
## Healthy indexes
HI.IN.PATH = snakemake@input[[4]]#'Classification/results/healthy_rf_results_all.RDS'
## PC1 scores
PCS.IN.PATH = snakemake@input[[5]]#"Integration_output/jive/subject/prcomp_list.rds"
## Microarray modules
MICROARRAY.MODULES.IN.PATH = snakemake@input[[6]]#'Data/Microarray/analysis_output/WGCNA/modules.rds'
## Microarray classifier gene signature
MICROARRAY.CLASSIFIER.IN.PATH = snakemake@input[[7]]#'Classification/transcriptional_surrogates/microarray_classifier_signatures.RDS'

## Transcriptional surrogate signatures for HI, PC1, and most significant features from classifier
SIGNATURES.OUT.PATH = snakemake@output[[1]]#'Classification/transcriptional_surrogates/surrogate_signatures.RDS'

# Load utility function
source('scripts/util/Signatures/create_signature.R')

# Initiate signatures list
all.signatures = list()

# Load data and metadata
matrices = readRDS(MATRICES.IN.PATH)
meta = readRDS(META.IN.PATH)

## Get the results of the pvalue permutations
results = readRDS(PVALS.IN.PATH)

# Extract the data used in the classifier, and all of the microarray features
X = matrices$all.modules.plus.grey.with.tbnks
Y = matrices$microarray.features

# Restrict data to only patients
#X = X[meta$condition != 'Healthy', ]
#Y = Y[meta$condition != 'Healthy', ]

# Remove the 'microarray.features.' prefix from the microarray feature names
colnames(Y) = gsub('microarray\\.features\\.', '', colnames(Y))

# First, we get the surrogate signatures for the healthy index driving features

## Get the feature pvalues for the full classifier 
result = results$all.modules.plus.grey.with.tbnks

## Get the features with an FDR of less than .2
features = names(result)[p.adjust(result, 'fdr') < .2]

## For each feature, make a microarray surrogate signature for that feature
signatures = lapply(features, function(feature) {util.make.signature(X[,feature], Y)})

## Name each suggoate signature based on the feature being approximated 
names(signatures) = features

## Add these signatures to the list of signatures 
all.signatures = append(all.signatures, signatures)

# Next, we make a surrogate signature for the healthy index itself

## Read in the healthy indexes
healthy.indexes = readRDS(HI.IN.PATH)

## Subset to only patients
#healthy.indexes = healthy.indexes[meta$condition != "Healthy", ]

## Choose the healthy index from the full classifier
healthy.index = healthy.indexes$all.modules.plus.grey.with.tbnks

## Make the healthy index names the patient ids
names(healthy.index) = rownames(healthy.indexes)

## Ensure that the healthy index subjects are in the same order as the microarray features
stopifnot(all(names(healthy.index) == rownames(Y)))

## Construct the healthy index singature
signature = util.make.signature(healthy.index, Y)

## Add the healthy index signature to the list of signatures
all.signatures[['healthy.index']] = signature

# Now we make a surrogate signature for the PC1 Signature
## Extract the jPC1 scores
jive = readRDS(PCS.IN.PATH)
joint = jive$joint$x
jPC1 = joint[,1]

## Subset and reorder the jPC1 scores to correspond to the subjects / order 
## in the microarray features
jPC1 = jPC1[rownames(Y)]

## Construct the jPC1 surrogate signature
signature = util.make.signature(jPC1, Y)

## Add this signature to the list of signatures
all.signatures[['PC1']] = signature

# Now we add in pseudo-signatures of the microarray modules. As these do not
# require data type conversions, we simply let them be composed of all the genes
# in a given module

## We extract the microarray modules that with significant GVIs in the microarray
## module classifier
result = results$microarray.modules
modules = names(result)[p.adjust(result, 'fdr') < .2]

## We load the microarray modules memberships
module.memberships = readRDS(MICROARRAY.MODULES.IN.PATH)

## For each microarray module that passed the significance threshold
signatures = lapply(modules, function(module) {
  ## We change the name of the module's constituent genes
  ## to eliminate the microarray.modules prefix
  module = gsub('microarray\\.modules\\.', '', module)
  
  ## We get the genes belonging to the module
  module.members = names(module.memberships)[module.memberships == module]
  
  ## We return the pseudo-signature
  list(positive = module.members, negative = NULL)
})

## We name the signatures accordingly
names(signatures) = modules

## And add them to our list of signatures
all.signatures = append(all.signatures, signatures)

# Finally, we load and add the microarray logistic regression classifier top features
signature = readRDS(MICROARRAY.CLASSIFIER.IN.PATH)
all.signatures[['microarray.classifier']] = signature

# We remove any empty signatures we may have picked up
all.signatures = all.signatures[sapply(all.signatures, function(signature) {length(unlist(signature))}) > 0]

## And we save the signatures
saveRDS(all.signatures, SIGNATURES.OUT.PATH)
