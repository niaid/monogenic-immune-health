require(Biobase)

get_module_scores = function(eset.train, eset.test, modules) {
  
  # extract the data from the esets
  X.train = t(exprs(eset.train))
  X.test = t(exprs(eset.test))
  
  # get module scores
  scores = mapply(function(module, module_name) {
    # Subset testing and training data to the features in the module
    X.train.module = X.train[, module, drop = FALSE]
    X.test.module = X.test[, module, drop = FALSE]
    
    # Train pca on the subsetted training data
    pca = prcomp(X.train.module, scale. = TRUE, center = TRUE, rank. = 1)
    direction = ifelse(mean(sign(cor(pca$x, X.train.module))) > 0, 1, -1)
    variance.explained = summary(pca)$importance[[2,1]]
    print(paste0('Variance explained by ', module_name,' module: ', as.character(variance.explained)))
    
    # Test pca on the subsetted training data
    module.scores = predict(pca, X.test.module) * direction
  }, modules, names(modules))
  
  # Name the rows and columns
  colnames(scores) = names(modules)
  rownames(scores) = rownames(X.test)
  
  # Make trivial feature meta data
  feature.meta = data.frame(module_name = names(modules))
  rownames(feature.meta) = names(modules)
  
  # Wrap the scores up into its own eset
  scores = ExpressionSet(t(scores))
  phenoData(scores) = phenoData(eset.test)
  featureData(scores) = AnnotatedDataFrame(feature.meta)
  
  return(scores)
}