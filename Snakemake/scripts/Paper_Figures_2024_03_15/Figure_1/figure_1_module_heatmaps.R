# This script creates the feature-feature correlation heatmaps with annotations for the gene and protein modules.
# Feature-feature correlations are computed using subject-level data and just plotted for the stable features.
# The subsetting is done primarily in order to save time when plotting the microarray heatmap.

# Load libraries
library(ggplot2)
library(ggrepel)
library(reshape2)
library(Biobase)
library(pheatmap)

# Load utility functons
source('scripts/util/Plotting/feature_and_module_heatmaps.R')
source('scripts/util/Processing/averageRepeatSamples.R')

# Set paths
SOMALOGIC.ESET.IN.PATH = snakemake@input[[1]]#'Data/Somalogic/data_analysis_ready/analysis_ready_subject_level_training_somalogic.rds'
SOMALOGIC.VP.IN.PATH = snakemake@input[[2]]#'Data/Somalogic/analysis_output/stability/somalogic_features_standard_vp.rds'
SOMALOGIC.MODULES.IN.PATH = snakemake@input[[3]]#'Data/Somalogic/analysis_output/wgcna_results/modules.rds'
MICROARRAY.ESET.IN.PATH = snakemake@input[[4]]#'Data/Microarray/data_analysis_ready/eset_batch_training_subject.rds'
MICROARRAY.VP.IN.PATH = snakemake@input[[5]]#'Data/Microarray/analysis_output/stability/microarray_features_standard_vp.rds'
MICROARRAY.MODULES.IN.PATH = snakemake@input[[6]]#'Data/Microarray/analysis_output/WGCNA/modules.rds' 

PROTEIN.HEATMAP.OUT.PATH = snakemake@output[[1]]#'Paper_1_Figures/Figure_1/1d_proteomic.png'
GENE.HEATMAP.OUT.PATH = snakemake@output[[2]]#'Paper_1_Figures/Figure_1/1d_transcriptomic.png'

# Create the feature-feature correlation heatmap for the proteins (only using stable features)
png(PROTEIN.HEATMAP.OUT.PATH, res = 72 * 3, height = 3*480, width = 3*480)
eset = readRDS(SOMALOGIC.ESET.IN.PATH)
vp = readRDS(SOMALOGIC.VP.IN.PATH)
modules = readRDS(SOMALOGIC.MODULES.IN.PATH)
features = rownames(vp)[vp$Residuals < .5]
eset = eset[features,]
modules = modules[features]
module.colors = setdiff(unique(modules), 'grey')
modules = factor(modules, levels = c(module.colors, 'grey'))
plot_feature_correlations(eset, modules)
dev.off()

# Create the feature-feature correlation heatmap for the genes (only using stable features)
png(snakemake@output[[2]], res = 72 * 3, height = 3*480, width = 3*480)
eset = readRDS(MICROARRAY.ESET.IN.PATH)
vp = readRDS(MICROARRAY.VP.IN.PATH)
modules = readRDS(MICROARRAY.MODULES.IN.PATH)
features = rownames(vp)[vp$Residuals < .5]
eset = eset[features,]
modules = modules[features]
module.colors = setdiff(unique(modules), 'grey')
modules = factor(modules, levels = c(module.colors, 'grey'))
plot_feature_correlations(eset, modules)
dev.off()

