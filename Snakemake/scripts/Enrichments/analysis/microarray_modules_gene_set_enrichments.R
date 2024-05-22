## Gets gene set enrichments for the microarray gene modules

# Source utility functions
source('scripts/util/Enrichment/hyperGeo.R')

# Set globals
## Microarray module memberships path
MODULES.IN.PATH = snakemake@input[[1]]#'Data/Microarray/analysis_output/WGCNA/modules.rds'
## List of gene sets paths
GENE.SETS.IN.PATH = snakemake@input[[2]]#'Gene_sets/processed/combined_gene_sets.RDS'

## Place to save microarray module enrichments
ENRICHMENTS.OUT.PATH = snakemake@output[[1]]#'Data/Microarray/analysis_output/enrichments/microarray_module_gene_set_enrichments.RDS'

# Load data
modules = readRDS(MODULES.IN.PATH)
gene.sets = readRDS(GENE.SETS.IN.PATH)

# Get set of all genes
universe = names(modules)

# Get all module colors
module_colors = unique(modules)

# Apply enrichment function to all modules
## For each module
enrichments = lapply(module_colors, function(module) {
  ## Make the hits the set of genes in the module
  hits = names(modules)[modules == module]
  ## Test to see if any gene set is enriched for these hits
  multiHyperGeoTests(gene.sets, universe, hits, minGeneSetSize = 5, pAdjustMethod = "BH")
})

# Name the enrichments based upon their corresponding module
names(enrichments) = module_colors

# Save results
saveRDS(enrichments, ENRICHMENTS.OUT.PATH)
