## Here we get the scores among the meta-analysis subjects
## using the surrogate signatures derived in the monogenic cohort.
## We do this both across diseases and in a single disease. Note that
## we only use the across disease scores, but in case reviewers/others
## are interested we also keep the disease-specific ones.

# Set seed
set.seed(1119)

# Load libraries and source utilities
library(MetaIntegrator)
source('scripts/util/Signatures/get_signature_scores.R')

# Set globals
## Cleaned cgps
CGPS.IN.PATH = snakemake@input[[1]]#'Reference/jamboree/data_analysis_ready/cgps_clean.RDS'
## Meta-analysis object
META.OBJ.IN.PATH = snakemake@input[[2]]#'Reference/jamboree/data_analysis_ready/meta_studies.RDS'
## Gene surrogate signatures for the features we wish to test
SIGNATURES.IN.PATH = snakemake@input[[3]]#'Classification/transcriptional_surrogates/surrogate_signatures.RDS'

## Meta-analysis study gene signature scores
SCORES.OUT.PATH = snakemake@output[[1]]#'Reference/jamboree/analysis_output/meta_integrator_signature_scores.RDS'
## Meta integrator gene signature scores for each disease subset
SUBSET.SCORES.OUT.PATH = snakemake@output[[2]]#'Reference/jamboree/analysis_output/meta_integrator_subset_signature_scores.RDS'

# Load data
cgps = readRDS(CGPS.IN.PATH)
metaObj = readRDS(META.OBJ.IN.PATH)
signatures = readRDS(SIGNATURES.IN.PATH)

# Instantiate a function to turn a meta object into a meta score object
get_scores = function(metaObj, signatures) {
  ## Extract the data objects from the metaObj
  dataObjs = metaObj$originalData
  ## For each set of data object
  scores = lapply(dataObjs, function(dataObj) {
    ## Get the signature scores for each signature and turn it into a
    ## new expression matrix
    expr = t(dataObj$expr)
    expr = util.get_signature_scores(expr, signatures) 
    ## Replace the data object with the new expression matrix
    dataObj$expr = t(expr)
    ## Replace the names of the original features with the names of the signatures
    dataObj$keys = rownames(dataObj$expr)
    ## Check that this data object has the correct format
    stopifnot(checkDataObject(dataObj, "Dataset"))
    return(dataObj)
  })
  ## Wrap the scores into a meta object
  metaObj = list(originalData = scores)
  ## Check that this meta object has the correct format
  stopifnot(checkDataObject(metaObj, "Meta", "Pre-Analysis"))
  return(metaObj)
}

# We get the signature scores across diseases
metaObj.scores = get_scores(metaObj, signatures)

# We get the signature scores in each disease
diseases = names(cgps)
## For each disease
metaObj.subset.scores = lapply(diseases, function(disease) {
  ## Find the studies for a single disease
  studies = names(cgps[[disease]])
  ## Subset the meta object to these studies
  metaObj.subset = list(originalData = metaObj$originalData[studies])
  ## Get the signature scores for these studies
  scores = get_scores(metaObj.subset, signatures)
})

names(metaObj.subset.scores) = diseases

# We save the results
saveRDS(metaObj.scores, SCORES.OUT.PATH)
saveRDS(metaObj.subset.scores, SUBSET.SCORES.OUT.PATH)