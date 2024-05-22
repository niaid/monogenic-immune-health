require(WGCNA)
require(pheatmap)
require(RColorBrewer)

plot_feature_correlations = function(es, colors, method = 'pearson', main = 'Feature Correlation Heatmap') {
  
  ## Get the correlation matrix from the expression set
  cor.mat = cor(t(exprs(es)), method = method)
  
  ## Ensure that the colors are a factor
  #colors = factor(colors)
  
  ## Reorder the features so that they are together
  feature.order = order(colors)
  cor.mat = cor.mat[feature.order, feature.order]
  colors = colors[feature.order]
  
  ## Create the annotations for the rows and columns
  annotation_col = data.frame(module = colors)
  rownames(annotation_col) = colnames(cor.mat)
  annotation_row = data.frame(module = colors)
  rownames(annotation_row) = rownames(cor.mat)
  
  ## Give colors to each module
  module.colors = unique(colors)
  names(module.colors) = module.colors
  annotation_colors = list(module = module.colors)
  
  breaksList = seq(-1, 1, by = .01)
  
  p = pheatmap(cor.mat, annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = annotation_colors,
               annotation_legend = FALSE, show_rownames = FALSE, show_colnames = FALSE, cluster_rows = FALSE, cluster_cols = FALSE, 
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList,
               main = main)
  print(p)
}

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

plot_sample_module_scores = function(es, scores, palette, main = 'Sample Module Scores (Normalized)') {
  
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

plot_feature_value_and_module_score_correlations = function(es, scores, main = 'Feature Value and Module Score Correlations') {
  module_cor = cor(t(exprs(es)), scores)
  module_colors = colnames(module_cor)
  module_colors = factor(module_colors)
  names(module_colors) = module_colors
  annotation_col = data.frame(module = module_colors)
  rownames(annotation_col) = colnames(module_cor)
  annotation_colors = list(module = module_colors)
  
  p = pheatmap(module_cor, annotation_col = annotation_col, annotation_colors = annotation_colors, annotation_legend = FALSE,
               show_rownames = FALSE, show_colnames = FALSE, cluster_cols = FALSE, main = main)
  print(p)
}