# Plots heatmaps 

# Load libraries
library(limma)
library(Biobase)
library(dplyr)
library(ComplexHeatmap)

# Set paths
SEX.LINKED.DE.RESULT.IN.PATHS = list(
  somalogic.modules = snakemake@input[[1]],#'Data/Somalogic/analysis_output/sex_related_de_signatures/somalogic_modules_DE_results.rds',
  somalogic.features = snakemake@input[[2]],#'Data/Somalogic/analysis_output/sex_related_de_signatures/somalogic_features_DE_results.rds',
  microarray.modules = snakemake@input[[3]],#'Data/Microarray/analysis_output/sex_related_de_signatures/microarray_modules_DE_results.rds',
  microarray.features = snakemake@input[[4]],#'Data/Microarray/analysis_output/sex_related_de_signatures/microarray_features_DE_results.rds',
  tbnks.features = snakemake@input[[5]]#'Data/TBNK/analysis_output/sex_related_de_signatures/tbnks_features_DE_results.rds'
)

FIGURE.OUT.PATH = snakemake@output[[1]]#'Paper_1_Figures/Supplemental_Figure_2/sex_condition_interaction_effects.pdf'

# Load data
resultss = lapply(SEX.LINKED.DE.RESULT.IN.PATHS, readRDS)

# Iniate the pdf
pdf(FIGURE.OUT.PATH)

# For each data type 
for(data.type in names(resultss)) {
  
  # Get the DE sex-condition interaction results for that data type
  results = resultss[[data.type]]
  
  # Get the adjusted pvalues associated with each condition for the interaction term
  result = sapply(results, function(result) {result$adj.P.Val[,ncol(result$adj.P.Val)]})
  
  # Get the effect sizes associated with each condition for the interaction term
  effects = sapply(results, function(result) {result$effect.size[,ncol(result$effect.size)]})
  
  # Change pvalues to signed log10 pvalues (signed based upon direction of the effect in female patients compared to healthy males
  # after correcting for condition and gender separately)
  result = -log10(result)
  
  # Attach a sign to the pvalue corresponding to the sign of the 
  result = result * sign(effects)
  
  # If there are more than 30 features
  if(nrow(result) > 30) {
    # Subset to only features that are significant at a .05 adjusted pvalue cutoff
    result = result[rowMax(result) > -log10(.05) | rowMin(result) < log10(.05), , drop = F]
    data.type = paste0('Negative log10\nsigned adjusted pvalues:\n',
                       data.type,
                       '\n(only significant\nfeatures shown;\npositive indicates\nhigher compared\nto healthy males\nafter correcting for \ncondition and gender\nseparately.)')
  } else{
    data.type = paste0('Negative log10\nsigned adjusted pvalues:\n',
                       data.type,
                       '\n(positive indicates\nhigher compared\nto healthy males\nafter correcting for \ncondition and gender\nseparately.)')
  }
  
  # Make the heatmap
  h = Heatmap(result, name = data.type, row_names_gp = gpar(fontsize = min(450/nrow(result), 12)))
  
  # And print it
  print(h)
}
dev.off()