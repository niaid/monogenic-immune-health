# Make figure 2 papers without patient who have been on IFN gamma

# Set seed
set.seed(150)

# Load libraries and source utilities function
library(ComplexHeatmap)
library(Biobase)
library(circlize)
source('scripts/util/Plotting/DE_heatmap.R')

# Set paths
ESETS.IN.PATHS = list(
  somalogic.modules = snakemake@input[[1]],#'Data/Somalogic/analysis_output/wgcna_results/scores_subject_level.rds',
  somalogic.features = snakemake@input[[2]],#'Data/Somalogic/data_analysis_ready/analysis_ready_subject_level_training_somalogic.rds',
  microarray.modules = snakemake@input[[3]],#'Data/Microarray/analysis_output/WGCNA/array_subject_scores.rds',
  microarray.features = snakemake@input[[4]],#'Data/Microarray/data_analysis_ready/eset_batch_training_subject.rds',
  tbnks = snakemake@input[[5]]#'Data/TBNK/data_analysis_ready/tbnk_training_subject_level_eset.rds'
)

RESULTS.IN.PATHS = list(
  somalogic.modules = snakemake@input[[6]],#'Data/Somalogic/analysis_output/differential_expression/somalogic_modules_DE_results_no_gamma.rds',
  somalogic.features = snakemake@input[[7]],#'Data/Somalogic/analysis_output/differential_expression/somalogic_features_DE_results_no_gamma.rds',
  microarray.modules = snakemake@input[[8]],#'Data/Microarray/analysis_output/differential_expression/microarray_modules_DE_results_no_gamma.rds',
  microarray.features = snakemake@input[[9]],#'Data/Microarray/analysis_output/differential_expression/microarray_features_DE_results_no_gamma.rds',
  tbnks = snakemake@input[[10]]#'Data/TBNK/analysis_output/differential_expression/tbnks_features_DE_results_no_gamma.rds'
)

PROTEIN.MODULES.IN.PATH = snakemake@input[[11]]#'Data/Somalogic/analysis_output/wgcna_results/modules.rds'
GENE.MODULES.IN.PATH = snakemake@input[[12]]#'Data/Microarray/analysis_output/WGCNA/modules.rds'
MEDICATIONS.IN.PATH = snakemake@input[[13]]#'Medications/medications.types.rds'

MAIN.HEATMAPS.OUT.PATH = snakemake@output[[1]]#"Paper_1_Figures/Figure_2/feature_heatmaps.pdf"
GENE.MODULE.HEATMAPS.OUT.PATH = snakemake@output[[2]]#"Paper_1_Figures/Figure_2/gene_heatmaps_by_module.pdf"
PROTEIN.MODULE.HEATMAPS.OUT.PATH = snakemake@output[[3]]#"Paper_1_Figures/Figure_2/protein_heatmaps_by_module.pdf"

# Load data
esets = lapply(ESETS.IN.PATHS, readRDS) # Load expression set data
results = lapply(RESULTS.IN.PATHS, readRDS) # Load results from limma DE testing
meds = readRDS(MEDICATIONS.IN.PATH) # Load medications

# Remove patients who were on gamma at any time point
## Get the visits during which patients were on subject
gamma.subjects = unique(meds$patient_id[meds$IFN.gamma])

## Remove the subjects from all the esets
esets = lapply(esets, function(eset) { eset[, !eset$patient_id %in% gamma.subjects] })

# Get the patient-condition pairings from each eset in order to calculate the number of subjects coming from each data type
sample.datas = lapply(esets, function(eset) {
  sample.data = pData(eset)
  sample.data[, c('patient_id','condition')]
})

# Get the condition counts associated with each data type
condition_counts = lapply(sample.datas, function(sample.data) {
  counts = table(sample.data$condition)
  counts = counts[names(counts) != "Healthy"]
})

# Make settings for all the plots

## This list controls what the heatmaps will be labeled as
result_names = list(
  somalogic.modules = 'Protein Modules (no gamma)\n',
  somalogic.features = 'Protein Features (no gamma)\n',
  microarray.modules = 'Gene Modules (no gamma)\n',
  microarray.features = 'Gene Features (no gamma)\n',
  tbnks = 'TBNK Features (no gamma)\n',
  all = 'Modules + TBNKs (no gamma)\n'
)

## This list controls the sizes of the *'s associated with significant features in the heatmap boxes
pnt_szs = list(
  somalogic.modules = 20,
  somalogic.features = 0,
  microarray.modules = 20,
  microarray.featues = 0,
  tbnks = 7,
  all = 3
)

## This list controls the sizes of the labels for the heatmap rows and columns
fnt_szs = list(
  somalogic.modules = 12,
  somalogic.features = 10,
  microarray.modules = 12,
  microarray.featues = 12,
  tbnks = 12,
  all = 8
)

# Plot the feature by feature heatmaps

## Open a pdf document
pdf(MAIN.HEATMAPS.OUT.PATH, width = 10, height = 10)

## For each data type
mapply(function(result, name, pnt_sz, fnt_sz, condition_count) {
  ## Get the associated statistics for the versus all comparison
  stats.all = result$versus.all
  
  ## Get the associated statistics for the versus healthy comparison
  stats.healthy = result$versus.healthy
  
  ## Plot the heatmap for the t-statistics associated with the DE signatures for that data type 
  util.plot.de.heatmap(stats.healthy, stats.all, name, condition.count = condition_count,
                       pnt.sz = pnt_sz, fnt.sz = fnt_sz, stat = 't')
}, results, result_names, pnt_szs, fnt_szs, condition_counts)
## Close the pdf
dev.off()

# Plot the DE results associated with the microarray features within each gene module

## Open a PDF
pdf(GENE.MODULE.HEATMAPS.OUT.PATH, width = 10, height = 10)
## Get the results of the DE testing for microarray features
result = results$microarray.features
## Get the microarray module associated with each feature
modules = readRDS(GENE.MODULES.IN.PATH)
## Get the vector of all microarray modules
colors = unique(modules)
## For each of the microarray modules
mapply(function(color) {
  print(color)
  ## Get the features associated with that module
  module = names(modules)[modules == color]
  ## Get the names of all features from the DE testing
  features = rownames(result$versus.all$effect.size)
  ## Subset these features to those in the module
  features = features[features %in% module]
  ## If there are no features in the module
  if(length(features) < 1) {
    ## Continue without plotting
    print("nothing to do")
    return(0)
  }
  ## Otherwise, proceed to the plotting
  ## Get the statistics associated with the versus all comparison
  stats.all = lapply(result$versus.all, function(x) {x[features, , drop = FALSE]})
  ## Get the statistics associated with the versus healthy conparison
  stats.healthy = lapply(result$versus.healthy, function(x) {x[features, , drop = FALSE]})
  ## Name the heatmap based on the color of the module
  name = paste0(color, ' gene module')
  ## Set the '*' size
  pnt_sz = min(round(120 / length(features)), 20)
  ## Set the label font size
  fnt_sz = min(ceiling(320 / length(features)), 12)
  ## Plot the heatmap
  util.plot.de.heatmap(stats.healthy, stats.all, name, condition.count = condition_counts$microarray.features,
                       pnt.sz = pnt_sz, fnt.sz = fnt_sz, stat = 't')
}, colors)
## Close the PDF
dev.off()

# Plot the DE results associated with the microarray features within each gene module

## Open a PDF
pdf(PROTEIN.MODULE.HEATMAPS.OUT.PATH, width = 10, height = 10)
## Get the results of the DE testing for somalogic features
result = results$somalogic.features
## Get the somalogic module associated with each feature
modules = readRDS(PROTEIN.MODULES.IN.PATH)
## Subset these features to those in the module
colors = unique(modules)
## For each of the somalogic modules
mapply(function(color) {
  print(color)
  ## Get the features associated with that module
  module = names(modules)[modules == color]
  ## Get the names of all features from the DE testing
  features = rownames(result$versus.all$effect.size)
  ## Subset these features to those in the module
  features = features[features %in% module]
  ## If there are no features in the module
  if(length(features) < 1) {
    ## Continue without plotting
    print("nothing to do")
    return(0)
  }
  ## Otherwise, proceed to the plotting
  ## Get the statistics associated with the versus all comparison
  stats.all = lapply(result$versus.all, function(x) {x[features, , drop = FALSE]})
  ## Get the statistics associated with the versus healthy conparison
  stats.healthy = lapply(result$versus.healthy, function(x) {x[features, , drop = FALSE]})
  ## Name the heatmap based on the color of the module
  name = paste0(color, ' protein module')
  ## Set the '*' size
  pnt_sz = min(round(120 / length(features)), 20)
  ## Set the label font size
  fnt_sz = min(ceiling(320 / length(features)), 12)
  ## Plot the heatmap
  util.plot.de.heatmap(stats.healthy, stats.all, name, condition.count = condition_counts$somalogic.features, 
                       pnt.sz = pnt_sz, fnt.sz = fnt_sz, stat = 't')
}, colors)
dev.off()

