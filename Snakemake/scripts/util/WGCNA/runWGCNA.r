library(WGCNA)
library(Biobase)

runWGCNA <- function(expressionset, out.dir = NA, method = 'hybrid', beta = 12, 
                     pamStage = TRUE, minModuleSize = 30, deepSplit = 1,
                     pamRespectsDendro = TRUE, diagnostic.plots.path = NULL, 
                     intermediate.results.path = NULL) {
  
  ## runs WGCNA algorithm and saves results in specified output directory
  
  ## expressionset is an expressionset object from the Biobase package
  ## outdir: the output directory where you would like the results and plots saved
  ## method: the method used for cutreeDynamic
  ## beta: the softpower. 12 is the suggested for a signed network; 
  ##      however, specify "auto" if you would like to automatically estimate this.
  ## pamStage: also used in cutreeDynamic
  ## minModuleSize: the minimum number of features that can be put into a module
  ## deepSplit: parameter used to control the tree cutting -- passed to WGCNA
  ## pamRespectsDendro: parameter used to control the PAM stage -- passed to WGCNA
  ## diagnostic.plots.path: where to save the output plots from WGCNA
  ## intermediate.results.path: where to save the intermediate objects from WGCNA
  
  ## The default parameters specified above are  essentially the default parameters
  ## for creating a signed network. 
  ## This parameter combination seems to work best for microarray,
  ## and other data types likely will need to try different parameter combinations
  
  ## intermediate results include the adjacency matrix and 1 - dissTOM. Set intermediate.results.path to NULL to avoid saving these
  ## Also can set diagnostic.plots.path to NULL to avoid saving diagnostic plots
  
  ## Set up
  
  datExpr <- t(exprs(expressionset))

  if(!is.null(diagnostic.plots.path)) {
    pdf(file = diagnostic.plots.path)
  }
  
  ## Soft thresholding
  
  # Set powers to investigate 
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed")
  
  if(!is.null(diagnostic.plots.path)) {
    # Plot the results:
    par(mfrow = c(1,2));
    cex1 = 0.9;
    
    # Scale-free topology fit index as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2", type="n", main = paste("Scale independence"));
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
    abline(h=0.90,col="red") # this line corresponds to using an R^2 cut-off of h
    
    # Mean connectivity as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  }
  
  ## Make adjacency and TOM
  
  if(beta == "auto") {
    softPower = sft$powerEstimate
  } else{
    softPower = beta
  }
  
  adjacency = adjacency(datExpr, power = softPower, type ="signed")
  
  TOM = TOMsimilarity(adjacency, TOMType = "signed");
  dissTOM = 1-TOM
  
  ## Make gene tree
  
  geneTree = hclust(as.dist(dissTOM), method = "average");
  
  ## Cut the gene tree into modules
  
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree,
                              distM = dissTOM,
                              pamStage = pamStage, 
                              deepSplit = deepSplit,
                              pamRespectsDendro = pamRespectsDendro,
                              minClusterSize = minModuleSize, method = method);
  
  # Convert numeric lables into colors
  dynamicColors = labels2colors(dynamicMods)
  
  if(!is.null(diagnostic.plots.path)) {
    
    # Plot the dendrogram and colors underneath
    plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        main = "Gene dendrogram and module colors")
  }
  
  ## Calculate eigengenes
  MEList = moduleEigengenes(datExpr, colors = dynamicColors)
  MEs = MEList$eigengenes
  
  ## Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs);
  
  ## Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average");
  
  MEDissThres = 0.25
  
  if(!is.null(diagnostic.plots.path)) {
    plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
    abline(h=MEDissThres, col = "red")
  }
  
  ## Call an automatic merging function
  merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  
  # The merged module colors
  mergedColors = merge$colors;
  
  # Plot the dendrogram again with merged modules
  if(!is.null(diagnostic.plots.path)) {
    plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                        c("Dynamic Tree Cut", "Merged dynamic"),
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
    dev.off()
  }
  
  ## Get the final modules
  modules = mergedColors
  names(modules) = featureNames(expressionset)
  
  if(!is.null(intermediate.results.path)) {
    ## Save module colors and labels for use in subsequent parts
    intermediate.results = list(adjacency, TOM, geneTree)
    names(intermediate.results) = c('adjacency', 'TOM', 'geneTree')
    saveRDS(intermediate.results, file = intermediate.results.path)
  }
  
  return(modules)
}
