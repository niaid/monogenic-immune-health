# Load library
library(dplyr)

# Set paths
VALIDATION.RESULTS.IN.PATH = snakemake@input[[1]]#'Reference/jamboree/analysis_output/results/jamboree_gene_level_results.RDS'
TABLE.OUT.PATH = snakemake@output[[1]]#'Paper_1_Figures/Figure_5_Tables/figure_5_meta_analysis_table.txt'

# Load data
scores = readRDS(VALIDATION.RESULTS.IN.PATH)

# Save to table
write.table(scores, file = TABLE.OUT.PATH, sep = "\t", row.names = F, col.names = T)
