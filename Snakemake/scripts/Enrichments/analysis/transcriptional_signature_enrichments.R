## Derive transcriptional surrogate signature gene set enrichments from Microarray correlates of features

# Load libraries
library(Biobase)
source('scripts/util/Enrichment/hyperGeo.R')

# Set globals
## The gene surrogate signatures for features of interest
SIGNATURES.IN.PATH = snakemake@input[[1]]#'Classification/transcriptional_surrogates/surrogate_signatures.RDS'
## The microarray subject-level training eset
ESET.IN.PATH = snakemake@input[[2]]#'Data/Microarray/data_analysis_ready/eset_batch_training_subject.rds'
## The gene sets
GENE.SETS.IN.PATH = snakemake@input[[3]]#'Gene_sets/processed/combined_gene_sets.RDS'

## The gene surrogate signature enrichments among the gene sets
ENRICHMENTS.OUT.PATH = snakemake@output[[1]]#'Classification/transcriptional_surrogates/surrogate_enrichments.RDS'

# Load data
signatures = readRDS(SIGNATURES.IN.PATH)
microarray = readRDS(ESET.IN.PATH)
gene.sets = readRDS(GENE.SETS.IN.PATH)

# Get gene universe
universe = rownames(microarray)

# We get enrichments from all genes
results = lapply(signatures, function(signature) {
  result = lapply(signature, function(half.signature) {
    print(length(universe))
    multiHyperGeoTests(gene.sets, universe, half.signature, minGeneSetSize = 5)
  })
  names(result) = names(signature)
  return(result)
})

names(results) = names(signatures)

saveRDS(results, ENRICHMENTS.OUT.PATH)
