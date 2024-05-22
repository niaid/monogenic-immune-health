## Determine if the genes in each signature tend to have significantly higher
## effect sizes compared with genes not in the signature

# Set seed
set.seed(179)

# Load libraries
library(Biobase)
library(MetaIntegrator)

# Set paths
## Meta-analysis metaintegrator study object
META.OBJECT.IN.PATH = snakemake@input[[1]]#'Reference/jamboree/data_analysis_ready/meta_studies.RDS'
## Gene surrogate signatures
SIGNATURES.IN.PATH = snakemake@input[[2]]#'Classification/transcriptional_surrogates/surrogate_signatures.RDS'

## Results of the enrichment analysis among meta-analysis gene effects
RESULTS.ENRICHMENTS.OUT.PATH = snakemake@output[[1]]#'Reference/jamboree/analysis_output/results/jamboree_enrichment_results.RDS'
## Gene-level meta-analysis results
META.ANALYSIS.RESULTS.OUT.PATH = snakemake@output[[2]]#'Reference/jamboree/analysis_output/results/jamboree_gene_level_results.RDS'

# Signature enrichments
## Get the object with the meta-analysis data
metaObject = readRDS(META.OBJECT.IN.PATH)
## Get the gene surrogate signatures we're interested in
signatures = readRDS(SIGNATURES.IN.PATH)

# Make sure the meta analysis object is in the correct format for meta integrator
stopifnot(checkDataObject(metaObject, "Meta", "Pre-Analysis"))

# Run the meta integrator meta analysis
results = runMetaAnalysis(metaObject)

# Get the statistics from the results
pools = results$metaAnalysis$pooledResults

# Initialize an empty list for the results for each signature
ress = list()

# For each signature
for(signature in names(signatures)) {
  
  ## Get all the genes in that signature
  genes = unname(unlist(signatures[[signature]]))
  
  ## Create a map between gene names and their meta analysis (absolute) effect sizes
  effects = abs(pools$effectSize)
  names(effects) = rownames(pools)
  
  ## See if the genes in the signature tend to have higher (absolute) effect sizes
  res = wilcox.test(effects[names(effects) %in% genes], effects[!names(effects) %in% genes], alternative = 'greater')
  
  ## Add the results to the results list
  ress[[signature]] = res
}

# Save the results
saveRDS(ress, RESULTS.ENRICHMENTS.OUT.PATH)
saveRDS(pools, META.ANALYSIS.RESULTS.OUT.PATH)

