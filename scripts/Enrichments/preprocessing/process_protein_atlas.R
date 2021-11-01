## Create a tissue gene set from the human protein atlas downloaded files

# Set paths
## File containing the human protein atlas information
ATLAS.IN.PATH = snakemake@input[[1]]#'Gene_sets/raw/Human_protein_atlas/proteinatlas.tsv'
## Gene sets for the tissues
GENE.SET.OUT.PATH = snakemake@output[[1]]#'Gene_sets/processed/tissue_gene_sets.RDS'

# Load data frame from protein atlas
df = read.table(ATLAS.IN.PATH, header = TRUE, sep = '\t', comment.char = '', quote = '', stringsAsFactors = FALSE)

# Instantiate set making function for a given enrichment level
make_set = function(enrichment.levels) {
  
  ## Subset to only enrichments of the desired enrichment.levels
  df.specific = df[df$RNA.tissue.category %in% enrichment.levels, ]
  
  ## Get the set of all tissues in the HPA
  tissues = df.specific$RNA.TS.TPM
  ## Separate the tissues that may be enriched
  tissues = strsplit(tissues, '\\;')
  ## Parse these tissues to remove associated numerical values
  tissues = sapply(tissues, function(tissue) {
    tissue = sapply(strsplit(tissue, '\\:'), function(x) {x[[1]]})
  })
  ## Only take the unique tissues
  tissues = unique(unlist(tissues))
  
  ## Get the genes associated with each tissue
  gene_sets = lapply(tissues, function(tissue) {
    tissue = paste0(tissue, ':')
    ### Get the set of all genes that are specific to a tissue
    genes = df.specific$Gene[grepl(tissue, df.specific$RNA.TS.TPM)]
    ### Get the unique set of all these genes
    genes = unique(genes)
  })
  
  names(gene_sets) = tissues
  
  return(gene_sets)
}

# Set the enrichment levels we are interested in making gene sets for
levels = list(
  strict = c("Tissue enriched"),
  medium = c("Tissue enriched", "Tissue enhanced"),
  general = c("Tissue enhanced", "Tissue enriched", "Group enriched")
)

# Get the tissue enrichments for each set of enrichment levels
tissue_sets = lapply(levels, make_set)

# Save results
saveRDS(tissue_sets, GENE.SET.OUT.PATH)
