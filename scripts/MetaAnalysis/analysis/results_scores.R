## Runs a meta-integrator meta-analysis on the
## signature scores of samples across the OMiCC
## jamboree data seed.

# Set seed
set.seed(79)

# Load libraries
library(Biobase)
library(MetaIntegrator)

# Set globals
## Meta integrator signature scores for our signatures of interest
META.SCORES.IN.PATH = snakemake@input[[1]]#'Reference/jamboree/analysis_output/meta_integrator_signature_scores.RDS'

## Jamboree meta-anlyses z-score-based meta result statistics
SCORES.RESULTS.OUT.PATH = snakemake@output[[1]]#'Reference/jamboree/analysis_output/results/jamboree_z_score_results.RDS'

# Signature significance results
## Get the signature scores meta-analysis object
metaObject = readRDS(META.SCORES.IN.PATH)

# Check that the metaObject is in the correct format
stopifnot(checkDataObject(metaObject, "Meta", "Pre-Analysis"))

# Run meta-analysis
results = runMetaAnalysis(metaObject)

# Subset to the desired result statistics
results = results$metaAnalysis

# Save results
saveRDS(results, SCORES.RESULTS.OUT.PATH)