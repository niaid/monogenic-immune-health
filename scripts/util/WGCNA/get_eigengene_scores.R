library(Biobase)

get_eigengene_scores = function(eset.train, eset.test, modules) {
  ## Creates a new set of axes based upon the first principal component loadings
  ## of each module through performing PCAs on the subset of eset.test features in
  ## a given module. Then the samples in eset.test are projected onto these axes to
  ## obtain the first principal component score of the eset.test samples for each module.
  ## Inputs:
  ## eset.train - The eset from which to learn the first principal component loadings for each module
  ## eset.test - The eset for which we wish to derive the module PC1 scores 
  ## modules - a named character vector, where the name is the feature name (i.e. "IgE"), which must be
  ## a subset of the rownames of eset.train and and eset.test, and the value is the module name (i.e. "blue", "green")
  ## Outputs:
  ## an eset containing the module scores for each sample in eset.test and the pData from eset.test
  
  # rearrange testing eset and training eset so that the modules are in the same order
  modules = modules[featureNames(eset.train)]
  eset.test = eset.test[featureNames(eset.train), ]
  
  # extract the data from the esets
  X.train = t(exprs(eset.train))
  X.test = t(exprs(eset.test))
  
  # remove grey module
  grey = modules == 'grey'
  X.train = X.train[, !grey]
  X.test = X.test[, !grey]
  modules = modules[!grey]
  
  # factorize the modules
  modules = unname(modules)
  modules = factor(modules)
  
  # get the PC scores in eset.test from the PC loading derived using eset.train
  scores = sapply(levels(modules), function(module) {
    
    # Subset testing and training data to the features in the module
    X.train.module = X.train[, modules == module, drop = FALSE]
    X.test.module = X.test[, modules == module, drop = FALSE]
    
    # Train pca on the subsetted training data
    pca = prcomp(X.train.module, scale. = TRUE, center = TRUE, rank. = 1)
    direction = ifelse(mean(sign(cor(pca$x, X.train.module))) > 0, 1, -1)
    
    # Get pca scores on the subsetted testing data
    module.scores = predict(pca, X.test.module) * direction
    
    return(module.scores)
    
  })
  
  # Name the rows and columns of eset.test scores
  colnames(scores) = levels(modules)
  rownames(scores) = rownames(X.test)
  
  # Make trivial feature meta data to put with the scores
  feature.meta = data.frame(module_name = levels(modules))
  rownames(feature.meta) = levels(modules)
  
  # Wrap the scores up into its own eset
  scores = ExpressionSet(t(scores))
  phenoData(scores) = phenoData(eset.test)
  featureData(scores) = AnnotatedDataFrame(feature.meta)
  
  return(scores)
}

get_eigengene_variance_explained = function(eset, modules) {
  
  # Make sure the modules are in the same order as the features
  modules = modules[featureNames(eset)]
  
  # extract the data from the eset
  X = t(exprs(eset))
  
  # remove grey module
  grey = modules == 'grey'
  X = X[, !grey]
  modules = modules[!grey]
  
  # factorize the modules
  modules = unname(modules)
  modules = factor(modules)
  
  # get variances explained by each module
  variances = sapply(levels(modules), function(module) {
    
    # Subset testing and training data to the features in the module
    X.module = X[, modules == module, drop = FALSE]
    
    # Train pca on the subsetted training data
    pca = prcomp(X.module, scale. = TRUE, center = TRUE, rank. = 1)
    
    # Get the variance of the first PC
    variance = summary(pca)$importance[["Proportion of Variance", "PC1"]]
  })
  
  return(variances)
}
