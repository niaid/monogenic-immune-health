## Subset the proteomic surrogate signatures for the HI and joint PC1
## to only the somamers that should be stable between serum and plasma.
## Additionally, in preparation for applying this signature to other data sets,
## add a version the signatures using somalogic IDs rather than protein names. 

# Load libraries
library(Biobase)

# Set paths
## The healthy index serum proteomic surrogate signature
HI.IN.PATH = snakemake@input[[1]]#'Classification/proteomic_surrogates/healthy.index.surrogates.RDS'
## The PC1 serum proteomic surrogate signature
PC1.IN.PATH = snakemake@input[[2]]#'Classification/proteomic_surrogates/PC1.score.surrogates.RDS'
## The somalogic training-level eset
ESET.IN.PATH = snakemake@input[[3]]#'Data/Somalogic/data_analysis_ready/analysis_ready_sample_level_training_somalogic.rds'
## The somalogic somamer table with Plasma and Serum dilutions
REF.TABLE.IN.PATH = snakemake@input[[4]]#'Data/Somalogic/raw/v1/somamer_table.txt'

## The healthy index plasma surrogate signature
HI.SOMAMER.OUT.PATH = snakemake@output[[1]]#'Classification/proteomic_surrogates/healthy.index.plasma.surrogates.RDS'
## The healthy index plasma surrogate signature with names as sequence ids rather than protein names
HI.ID.OUT.PATH = snakemake@output[[2]]#'Classification/proteomic_surrogates/healthy.index.plasma.surrogate.somaId.RDS'
## The PC1 index plasma surrogate signature
PC1.SOMAMER.OUT.PATH = snakemake@output[[3]]#'Classification/proteomic_surrogates/PC1.plasma.surrogates.RDS'
## The PC1 index plasma surrogate signature with names as sequence ids rather than protein names
PC1.ID.OUT.PATH = snakemake@output[[4]]#'Classification/proteomic_surrogates/PC1.plasma.surrogate.somaId.RDS'

# Load data
healthy.index.signature = readRDS(HI.IN.PATH)
PC1.score.signature = readRDS(PC1.IN.PATH)
eset = readRDS(ESET.IN.PATH)

# Read in the reference table with somamer information
ref = read.table(REF.TABLE.IN.PATH, sep = '\t', header = TRUE, comment.char = '', quote = c())

# Get the somamer feature meta data from the eset
feature.meta = fData(eset)

# Create a map from somamer ids to features names
conversion = feature.meta$SomaId
names(conversion) = rownames(feature.meta)

# Create a backward map
inv_conversion = names(conversion)
names(inv_conversion) = conversion

# Instantiate a function to subset a signature to the desired targets
process_signature = function(signature) {
  ## For each half signature (positive and negative), convert the proteins name to somamer ids
  signature.converted = lapply(signature, function(x) {unname(unlist(conversion[x]))})
  
  ## Make sure that there are no somamers in the signature not included in the reference table
  stopifnot(setdiff(unlist(signature.converted), ref$SomaId) == character(0))
  
  ## Subset the reference table to just those with the same dilutions in serum and plasma
  ref.subset = ref[ref$CommonDilution == 'No',]
  
  ## Subset each half signature (positive and negative) to just porteins with the same dilutions in
  ## serum and plasma
  signature.converted = lapply(signature.converted, function(x) {setdiff(x, ref.subset$SomaId)})
  
  ## Return both half signatures to the original protein names rather than somamer ids
  signature = sapply(signature.converted, function(x) {unname(inv_conversion[x])})
  
  return(list(targets = signature, somaIds = signature.converted))
}

# Convert the serum HI signature to a plasma HI signature
signature = process_signature(healthy.index.signature)

# Save the results
saveRDS(signature$targets, HI.SOMAMER.OUT.PATH)
saveRDS(signature$somaIds, HI.ID.OUT.PATH)

# Convert the serum PC1 signature to a plasma PC1 signature
signature = process_signature(PC1.score.signature)

# Save the results
saveRDS(signature$targets, PC1.SOMAMER.OUT.PATH)
saveRDS(signature$somaIds, PC1.ID.OUT.PATH)
