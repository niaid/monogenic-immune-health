## We clear our workspace
rm(list = ls())

## We load the relevant libraries and utility functions
library(Biobase)
source('util/WGCNA/runWGCNA.r')
source('util/WGCNA/get_eigengene_scores.R')
source('util/Plotting/feature_and_module_heatmaps.R')
source('util/Plotting/colors.R')

## We load the relevant data
results = readRDS('../Classification/results/healthy_rf_results_2.RDS')
matrices = readRDS('../Classification/random_forest_design_matrices.RDS')

## We extract the most important features based upon the classifier
significant.features = lapply(results, function(result) {
  pvals = result$pvals
  names(pvals)[p.adjust(pvals, method = 'fdr') <= .4]
})
significant.features = unique(unname(unlist(significant.features)))

## We extract (and subset the relevant matrices)
X = matrices$matrices$all.modules.plus.grey.with.tbnks
X = X[, significant.features]
Y = matrices$matrices$microarray.features
meta = matrices$meta

## We find all microarray features that are correlated at least .4 with the significant features
cor.mat = cor(Y, X)
max.cors = rowMax(cor.mat)
correlated.features = rownames(cor.mat)[max.cors >= .5]

## We make a matrix with just the correlated features
X.surrogates = Y[, correlated.features]

## We make a microarray surrogates expression set
feature.data = data.frame(features = colnames(X.surrogates))
rownames(feature.data) = colnames(X.surrogates)
subject.data = meta
eset = ExpressionSet(t(X.surrogates))
featureData(eset) = AnnotatedDataFrame(feature.data)
phenoData(eset) = AnnotatedDataFrame(subject.data)

## We get the modules and modules scores
modules = runWGCNA(eset, beta = 12)
scores = get_eigengene_scores(eset, eset, modules)

## Save the modules and scores
saveRDS(modules, '../Classification/transcriptional_surrogates/modules_.5.RDS')
saveRDS(eset, '../Classification/transcriptional_surrogates/signatures_.5.RDS')
saveRDS(scores, '../Classification/transcriptional_surrogates/scores_.5.RDS')