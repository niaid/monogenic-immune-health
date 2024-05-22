## Create a proteomic signature for the healthy index and joint PC1 using
## all of the stable proteomic features

# Load libraries and source utility functions
library(Biobase)
source('scripts/util/Signatures/create_signature.R')

# Set globals
## The healthy index for each subject
HI.IN.PATH = snakemake@input[[1]]#'Classification/results/healthy_rf_results_all.RDS'
## The jive PC1 for each subject
PC1.IN.PATH = snakemake@input[[2]]#'Integration_output/jive/subject/prcomp_list.RDS'
## The design matrices
MATRICES.IN.PATH = snakemake@input[[3]]#'Classification/design_matrices/healthy_all_design_matrices_all.RDS'

## The HI serum proteomic surrogate signature
HI.OUT.PATH = snakemake@output[[1]]#'Classification/proteomic_surrogates/healthy.index.surrogates.RDS'
## The PC1 serum proteomic surrogate signature
PC1.OUT.PATH = snakemake@output[[2]]#'Classification/proteomic_surrogates/PC1.score.surrogates.RDS'

# Load files
results = readRDS(HI.IN.PATH)
matrices = readRDS(MATRICES.IN.PATH)

# Get the somalogic features matrix
X = matrices$somalogic.features

# Remove the 'somalogic.features' prefix from the matrix column names
colnames(X) = gsub('somalogic\\.features\\.', '', colnames(X))

# Here, we get the signature for the HI:
# Get the healthy index for each subject 
preds = results$all.modules.plus.grey.with.tbnks
# Ensure that the rows of X correspond to the same subjects in the predictions vector
stopifnot(all(rownames(X) == rownames(preds)))
# Get the HI proteomic surrogate signature
HI.signature = util.make.signature(preds, X)

# Here, we do the same for the PC1 score:
# Extract the jPC1 score
jive = readRDS(PC1.IN.PATH)
joint = jive$joint$x
# Reorder the subjects so they correspond to the rows of the data matrix
joint = joint[rownames(X), , drop = FALSE]
# Choose only the PC1 score
joint = joint[, 'PC1']
# Get the PC1 proteomic surrogate signature
PC1.signature = util.make.signature(joint, X)

# Save results
saveRDS(HI.signature, HI.OUT.PATH)
saveRDS(PC1.signature, PC1.OUT.PATH)
