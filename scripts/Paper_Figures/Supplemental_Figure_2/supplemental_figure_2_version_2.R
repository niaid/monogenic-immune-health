# Load libraries
library(Biobase)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(dplyr)
library(reshape2)
library(pheatmap)

# Source utilities
source('scripts/util/Enrichment/proteinToGeneConversion.R')
source('scripts/util/Plotting/colors.R')

# Set paths
ESET.IN.PATH = snakemake@input[[1]]#'Data/Somalogic/analysis_output/stability/stable_somalogic_subject_level_features.rds'
COUNTS.IN.PATH = snakemake@input[[2]]#'Reference/protein_atlas/rna_tissue.tsv'

SOMALOGIC.MODULES.IN.PATH = snakemake@input[[3]]#'Data/Somalogic/analysis_output/wgcna_results/modules.rds'
MICROARRAY.MODULES.IN.PATH = snakemake@input[[4]]#'Data/Microarray/analysis_output/WGCNA/modules.rds'

SOMALOGIC.FEATURES.ESET.IN.PATH = snakemake@input[[5]]#'Data/Somalogic/data_analysis_ready/analysis_ready_sample_level_training_somalogic.rds'
MICROARRAY.FEATURES.ESET.IN.PATH = snakemake@input[[6]]#'Data/Microarray/data_analysis_ready/eset_batch_training_sample.rds'
TBNK.ESET.IN.PATH = snakemake@input[[7]]#'Data/TBNK/data_analysis_ready/tbnk_training_subject_level_eset.rds'
SOMALOGIC.MODULES.ESET.IN.PATH = snakemake@input[[8]]#'Data/Somalogic/analysis_output/wgcna_results/scores_subject_level.rds'
MICROARRAY.MODULES.ESET.IN.PATH = snakemake@input[[9]]#'Data/Microarray/analysis_output/WGCNA/array_subject_scores.rds'

SUPPLEMENTARY.FIGURE.2b.OUT.PATH = snakemake@output[[1]]#'Paper_1_Figures/Supplemental_Figure_2/S2b.pdf'
SUPPLEMENTARY.FIGURE.2c.OUT.PATH = snakemake@output[[2]]#'Paper_1_Figures/Supplemental_Figure_2/S2c.pdf'
SUPPLEMENTARY.FIGURE.2.EXTRA.OUT.PATH = snakemake@output[[3]]#'Paper_1_Figures/Supplemental_Figure_2/S2extra.pdf'

# Supplemental figure 2a -- breakdown of DE results by feature from protein grey module
## Find in Paper_1_Figures/Figure_2/protein_heatmaps_by_module.pdf

# Supplemental figure 2b/2c -- Normalized expression of grey module proteins / non-grey module proteins 
# in each tissue type from the Human Protein Atlas

## Load data
eset = readRDS(ESET.IN.PATH)
modules = readRDS(SOMALOGIC.MODULES.IN.PATH)
counts = read.table(COUNTS.IN.PATH, header = TRUE, sep = '\t', comment.char = '', quote = '')

## Get protein to gene map (eliminating multi-gene-mapping proteins)
protein_to_gene_map = make_protein_to_gene_map(eset)

## Get grey module proteins
grey.module = names(modules)[modules == 'grey']
grey.module = grey.module[grey.module %in% names(protein_to_gene_map)]

## Get coherent genes (those genes whose corresponding proteins are either all in or all not in the module)
coherent_genes = get_coherent_genes(protein_to_gene_map, grey.module)

## Convert the grey module proteins to genes 
grey.module = protein_to_gene_map[grey.module]
grey.module = intersect(grey.module, grey.module)
non.grey.module = setdiff(coherent_genes, grey.module)

## Subset tissue counts data to the coherent genes
counts = counts[counts$Gene.name %in% coherent_genes,]

## Make the counts matrix
counts = counts %>%
  dcast(Gene.name ~ Sample, value.var = 'Value') %>%
  as.data.frame()

rownames(counts) = counts$Gene.name

counts = counts %>%
  select(-Gene.name) %>%
  as.matrix()

## Scale the counts matrix across tissues
counts = counts %>% 
  t %>% 
  scale %>%
  t

## Split the genes into those corresponding to the grey module and those not corresponding to it
genes = rownames(counts)
group1 = genes %in% grey.module
group2 = genes %in% non.grey.module

## Restrict counts to group1 (genes in the grey module)
pdf(SUPPLEMENTARY.FIGURE.2b.OUT.PATH)
grey.counts = counts[group1, ]
grey.counts = grey.counts[complete.cases(grey.counts),]
h = Heatmap(
  grey.counts, show_row_names = F, 
  heatmap_legend_param = list(title = "Normalized\nExpression"),
  row_title = 'Grey Proteins'
)
draw(h)
dev.off()

## Restrict counts to group2
pdf(SUPPLEMENTARY.FIGURE.2c.OUT.PATH)
non.grey.counts = counts[group2, ]
non.grey.counts = non.grey.counts[complete.cases(non.grey.counts),]
h = Heatmap(
  non.grey.counts, show_row_names = F, 
  heatmap_legend_param = list(title = "Normalized\nExpression"),
  row_title = 'Non-Grey Proteins'
)
draw(h)
dev.off()

# Extra supplemental figures - sample level heatmaps of each data type
somalogic.module.colors = readRDS(SOMALOGIC.MODULES.IN.PATH)
microarray.module.colors = readRDS(MICROARRAY.MODULES.IN.PATH)

## Instantiate function to plot feature-level data
plot_sample_feature_values = function(es, colors, palette, main = 'Sample Feature Values (Normalized)') {
  
  ## Get the expression matrix from the expression set, with features as columns and samples as rows
  X = t(exprs(es))
  
  ## Scale the expression matrix
  X = scale(X)
  
  ## Get colors and conditions
  colors = factor(colors)
  conditions = factor(es$condition)
  
  ## Reorder the features so that they are together
  feature.order = order(colors)
  X = X[, feature.order]
  colors = colors[feature.order]
  
  ## Create the annotations for the rows and columns
  annotation_col = data.frame(module = colors)
  rownames(annotation_col) = colnames(X)
  annotation_row = data.frame(condition = conditions)
  rownames(annotation_row) = rownames(X)
  
  ## Create the list storing the colors associated with each annotation
  condition_colors = palette[1:length(levels(conditions))]
  names(condition_colors) = levels(conditions)
  module.colors = levels(colors)
  names(module.colors) = module.colors
  annotation_colors = list(module = module.colors, condition = condition_colors)
  
  ## Plot the heatmap
  p = pheatmap(X, annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = annotation_colors,
               annotation_legend = TRUE, show_rownames = FALSE, show_colnames = FALSE, cluster_cols = FALSE, main = main)
  print(p)
}

## Instantiate a function to plot module-level heatmaps
plot_sample_module_scores = function(es, palette, main = 'Sample Module Scores (Normalized)') {
  
  ## Extract the expression matrix
  scores = t(exprs(es))
  
  ## Scale the expression matrix
  scores = scale(scores)
  
  ## Get conditions
  rownames(scores) = colnames(es)
  conditions = factor(es$condition)
  
  ## Get colors
  module_colors = colnames(scores)
  names(module_colors) = module_colors
  
  ## Create the annotations for the rows and columns
  annotation_col = data.frame(module = module_colors)
  annotation_row = data.frame(condition = conditions)
  rownames(annotation_col) = colnames(scores)
  rownames(annotation_row) = rownames(scores)
  
  ## Create annotation colors list
  condition_colors = palette[1:length(levels(conditions))]
  names(condition_colors) = levels(conditions)
  annotation_colors = list(module = module_colors, condition = condition_colors)
  
  ## Plot the heatmap
  p = pheatmap(scores, annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = annotation_colors,
               annotation_legend = TRUE, show_rownames = FALSE, show_colnames = FALSE, cluster_cols = FALSE, main = main)
  
  print(p)
}

## Instantiate function for plotting TBNKs
plot_tbnks = function(es, palette, main = 'Sample Module Scores (Normalized)') {
  
  ## Extract the expression matrix
  scores = t(exprs(es))
  
  ## Scale the expression matrix
  scores = scale(scores)
  
  ## Get conditions
  rownames(scores) = colnames(es)
  conditions = factor(es$condition)
  
  ## Create the annotations for the rows and columns
  annotation_row = data.frame(condition = conditions)
  rownames(annotation_row) = rownames(scores)
  
  ## Create annotation colors list
  condition_colors = palette[1:length(levels(conditions))]
  names(condition_colors) = levels(conditions)
  annotation_colors = list(condition = condition_colors)
  
  ## Plot the heatmap
  p = pheatmap(scores, annotation_row = annotation_row, annotation_colors = annotation_colors,
               annotation_legend = TRUE, show_rownames = FALSE, show_colnames = TRUE, cluster_cols = FALSE, main = main)
  
  print(p)
}

somalogic = readRDS(SOMALOGIC.FEATURES.ESET.IN.PATH)
microarray = readRDS(MICROARRAY.FEATURES.ESET.IN.PATH)
tbnks = readRDS(TBNK.ESET.IN.PATH)
somalogic.modules = readRDS(SOMALOGIC.MODULES.ESET.IN.PATH)
microarray.modules = readRDS(MICROARRAY.MODULES.ESET.IN.PATH)

pdf(SUPPLEMENTARY.FIGURE.2.EXTRA.OUT.PATH, width = 8, height = 14)
plot_sample_feature_values(somalogic, somalogic.module.colors, high_contrast_palette())
plot_sample_feature_values(microarray, microarray.module.colors, high_contrast_palette())
plot_sample_module_scores(somalogic.modules, high_contrast_palette())
plot_sample_module_scores(microarray.modules, high_contrast_palette())
plot_tbnks(tbnks, high_contrast_palette())
dev.off()