## This script makes several gene sets to use for enrichments

# Source utility functions
source('scripts/util/Enrichment/gmtEnrichment.R')

# Set the paths
## For the GMTs to read in
GMT.IN.PATHS = list(
  reactome = snakemake@input[[1]],#'Gene_sets/raw/GMTs/reactome.gmt',
  go.bp = snakemake@input[[2]],#'Gene_sets/raw/GMTs/c5.bp.v6.2.symbols.gmt.txt',
  kegg = snakemake@input[[3]]#'Gene_sets/raw/GMTs/c2.cp.kegg.v6.2.symbols.gmt.txt'
)

## For the file with BTM information
BTMS.IN.PATH = snakemake@input[[4]]#'Gene_sets/raw/BTMs/btm_annotation_table.txt'

## For where to save the gene sets
GENE.SET.OUT.PATHS = list(
  reactome = snakemake@output[[1]],#'Gene_sets/processed/reactome.RDS',
  go.bp = snakemake@output[[2]],#'Gene_sets/processed/go.bp.RDS',
  kegg = snakemake@output[[3]],#'Gene_sets/processed/kegg.RDS',
  btms = snakemake@output[[4]]#'Gene_sets/processed/btm.rds'
) 

## For where to save a list with all the gene sets
GENE.SETS.OUT.PATH = snakemake@output[[5]]#'Gene_sets/processed/combined_gene_sets.RDS'

# Read in the GMTs
gene.sets = lapply(GMT.IN.PATHS, function(in.path) {
  gene.set = read.gmt(in.path)
})

# Get the btms

## Read in data
btm.dat <- read.table(BTMS.IN.PATH, header = TRUE, sep = '\t', comment.char = '', stringsAsFactors = FALSE)

## Save as list and give proper name
btms <- as.list(btm.dat$Module.member.genes)
btm.names <- paste(btm.dat$ID, btm.dat$Module.title, sep = "_")
names(btms) <- btm.names

## Remove commas and save as character vector where each gene is a component of the vector
btms <- lapply(btms, function(x) unlist(strsplit(x, split = ",")))

# Add the btms to the list of gene.sets
gene.sets[['btms']] = btms

# Save the gene sets separately
mapply(function(gene.set, out.path) {saveRDS(gene.set, out.path)}, gene.sets, GENE.SET.OUT.PATHS)

# Save the list of all gene sets
saveRDS(gene.sets, GENE.SETS.OUT.PATH)