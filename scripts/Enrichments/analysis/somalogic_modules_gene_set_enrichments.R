## Get gene set enrichments for somalogic proein modules. Because the gene sets are composed of genes
## and the proteomic modules are composed of proteins, we need to convert from gene to protein using the
## somamer information we have. However, sometimes, genes in a module are not 'coherent': sometimes
## a gene will correspond to multiple proteins, some of which fall into the module and some fall outside of it.
## When calculating the enrichments for the modules, we throw out the genes that are not coherent in that module.
## Proteins corresponding to multiple genes are taken out entirely.

# Load libraries
library(Biobase)

# Source utilities
source('scripts/util/Enrichment/hyperGeo.R')
source('scripts/util/Enrichment/proteinToGeneConversion.R')

# Set globals
## Somalogic module memberships
MODULES.IN.PATH = snakemake@input[[1]]#'Data/Somalogic/analysis_output/wgcna_results/modules.rds'
## Somalogic feature level eset, from which to get the fData for the somamers
ESET.IN.PATH = snakemake@input[[2]]#'Data/Somalogic/data_analysis_ready/analysis_ready_subject_level_training_somalogic.rds'
## The gene sets we wish to investigate for enrichments
GENE.SETS.IN.PATH = snakemake@input[[3]]#'Gene_sets/processed/combined_gene_sets.RDS'
## The tissue specific gene sets we wish to investigate for tissue enrichments
TISSUE.SETS.IN.PATH = snakemake@input[[4]]#'Gene_sets/processed/tissue_gene_sets.RDS'

## The location in which we wish to save the gene set enrichments for the somalogic
GENE.SETS.ENRICHMENTS.OUT.PATH = snakemake@output[[1]]#'Data/Somalogic/analysis_output/enrichments/somalogic_module_gene_set_enrichments.RDS'
## The location in which we wish to save the tissue enrichments for the somalogic
TISSUE.SETS.ENRICHMENTS.OUT.PATH = snakemake@output[[2]]#'Data/Somalogic/analysis_output/enrichments/somalogic_module_tissue_set_enrichments.RDS'

# Load data
modules = readRDS(MODULES.IN.PATH)
eset = readRDS(ESET.IN.PATH)
general.sets = readRDS(GENE.SETS.IN.PATH)
tissue.sets = readRDS(TISSUE.SETS.IN.PATH)

# Make a list combining tissue and gene sets for easy vectorization
gene.set.list = list(general = general.sets, tissue = tissue.sets)

# Get a map converting protein names to corresponding gene names based on the fData in somalogic
protein.to.gene.map = make_protein_to_gene_map(eset)

# Subset module to proteins in the map (i.e. those that correspond to only one gene)
modules = modules[names(modules) %in% names(protein.to.gene.map)]

# Get all the module colors
module_colors = unique(modules)

# Get the enrichments
## For both tissue gene sets and the normal gene sets
enrichments = lapply(gene.set.list, function(gene.sets) {
  ## For each module
  enrichments = lapply(module_colors, function(module_color){
    ## Get the proteins in the module
    module = names(modules)[modules == module_color]
    ## Get each gene for which all its corresponding proteins fall into or out of the module
    universe = get_coherent_genes(protein.to.gene.map, module)
    ## Extract the genes for which all corresponding proteins fall into the module
    hits = intersect(universe, protein.to.gene.map[module])
    ## Run gene set tests on the 'coherent' genes found above
    gene.sets.enrichments = multiHyperGeoTests(gene.sets, universe, hits, minGeneSetSize = 5, pAdjustMethod = "BH")
  })
  ## Name the enrichments based on the module
  names(enrichments) = module_colors
  ## Return result
  return(enrichments)
})

# Split enrichments into the general gene sets and the tissue gene sets and save each
saveRDS(enrichments$general, GENE.SETS.ENRICHMENTS.OUT.PATH)
saveRDS(enrichments$tissue, TISSUE.SETS.ENRICHMENTS.OUT.PATH)
