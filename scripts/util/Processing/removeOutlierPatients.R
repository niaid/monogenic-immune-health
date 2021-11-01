library(WGCNA)

removeOutlierPatients = function(eset, cutHeight, minSize = 10) {
  ## Uses hierarchical agglomerative clustering with a euclidean metric
  ## to remove outlier samples.
  ## Inputs:
  ## eset - the eset to analyze for outliers
  ## cutHeight - an integer specifying the height at which to cut the tree
  ## minSize - an integer specifying the number of samples specifying the minumum 
  ##           number of samples required to be in a cut branch for these samples
  ##           to be not considered outliers
  ## Outputs:
  ## An eset subsetted to the samples not considered to be outliers
  
  ## Extract the data 
  X = t(exprs(eset))
  
  ## Scale the data
  X = scale(X)
  
  ## Cluster samples
  sampleTree = hclust(dist(X), method = 'average')
  
  ## Plot clustering
  par(cex = .1, mar=c(5,10,5,1))
  plot(sampleTree, ylab = '', xlab = '', main = '', yaxt = "n")
  par(cex = 1)
  title(xlab = '', ylab = 'height', main = 'Dendrogram of Samples')
  axis(2, cex.axis = .5)
  
  ## Remove unwanted samples
  clust = cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = minSize)
  eset = eset[, clust == 1]
  X = X[clust == 1, ]
  
  ## Cluster samples again
  sampleTree = hclust(dist(X), method = 'average')
  
  ## Plot clustering again
  par(cex = .1, mar=c(5,10,5,1))
  plot(sampleTree, ylab = '', xlab = '', main = '', yaxt = "n")
  par(cex = 1)
  title(xlab = '', ylab = 'height', main = 'Dendrogam of Filtered Samples')
  axis(2, cex.axis = .5)
  
  ## Return reduced expression set
  return(eset)
}