hyperGeoTest = function(gene_set, universe, hits, fix_universe = TRUE, fix_gene_set = TRUE, fix_hits = TRUE) {
  ## Performs a single hyper geometric test using one gene_set and one set of hits
  ## Inputs:
  ## gene_set - a character vector with the names of all genes belonging to some set of interest
  ## universe - a character vector with the names of all genes to consider. Genes in 'gene_set' and in
  ##            'hits' will be subsetted to just those belonging to the universe
  ## hits - a character vector with the names of all genes we belonging to another set of interest
  ## fix_universe - it is required that each gene in the universe is only repeated once. If performing
  ##                this operation many times with the same universe, however, it may be inefficient to
  ##                ensure that each gene in the univese only appears once. Thus if fix_universe is TRUE,
  ##                the genes in the universe will be subset to the unique genes. 
  ##                If fix_universe is FALSE, it is trusted that the user has already ensured that each 
  ##                gene only appears once in 'universe'.
  ## fix_gene_set - similar to fix universe, it is required that each gene in the gene set is unique
  ##                and that all the genes in the gene_set belong to the universe. If fix_gene_set is TRUE,
  ##                the gene_set will be modified so that each gene only appears once, and all genes in the
  ##                gene set are part of the universe.
  ## fix_hits - analogous to fix_gene_set but for the hits
  ## Outputs:
  ## list with structure
  ## $p.value - numeric indicating the p-value associated with the hypergeometric test that genes
  ##             in the gene set are more commonly part of the hits than would be expected at random 
  ## $observed.hits - numeric indicating the number of genes in the gene set that belong to hits
  ## $expected.hits - numeric indicating the number of genes in the gene set that we would expect belong to hits
  ##                  under the null hypothesis of the hypergeometric test
  ## $odds.ratio - numeric indicating the odds that a gene belonging to the gene set would also be a hit
  ##               divided by the odds that a gene not belonging to the gene set would also be a hit
  
  # If the universe has not already been subset to unique genes
  if(fix_universe) {
    # Do so
    universe = unique(universe)
  }
  
  # If the hits are not already unique and a subset of the universe
  if(fix_hits) {
    # Make them so
    hits = unique(hits)
    hits = intersect(hits, universe)
  }
  
  # If the genes in the gene set are not already unique and a subset of the universe
  if(fix_gene_set) {
    # Make them so
    gene_set = unique(gene_set)
    gene_set = intersect(gene_set, universe)
  }
  
  # Get the number of genes in the genes set that are also hits
  x11 = length(intersect(hits, gene_set))
  # Get the number of hits that are not in the gene set
  x12 = length(hits) - x11
  # Get the number of genes in the gene set that are not hits
  x21 = length(gene_set) - x11
  # Get the number of genes in the universe that are neither hits nor in the gene set
  x22 = length(universe) - x11 - x12 - x21
  
  # Build a 2x2 contingency table with these counts
  X = matrix(c(x11, x12, x21, x22), ncol = 2, nrow = 2)
  
  # Test the hypothesis that genes in the gene set are also hits more
  # often than random
  result = fisher.test(X, alternative = 'greater')
  
  # Extract the pvalue associated with the result
  p.value = result$p.value
  
  # Get the number of hits observed in the gene set
  observed.hits = x11
  
  # Get the expected number of hits in the gene set under the null hypothesis of the hypergeometric test
  expected.hits = length(hits) / length(universe) * length(gene_set)
  
  # Get the odds ratio
  odds.ratio = (x11 * x22) / (x12 * x21)
  
  # Put together all these statsitics
  stats = list(
    'p.value' = p.value,
    'observed.hits' = observed.hits,
    'expected.hits' = expected.hits,
    'odds.ratio' = odds.ratio
  )
  
  return(stats)
}


## This function was an implementation of the HTSanalyzeR package's function,
## as their function had problems regarding installation issues on our singularity
## container. Please find the citation for this package below:
## Xin Wang and Camille Terfve and John C. Rose and Florian
## Markowetz. HTSanalyzeR: an R/Bioconductor package for integrated
## network analysis of high-throughput screens. Bioinformatics 27:6
## 879 (2011).
multiHyperGeoTest = function(collectionOfGeneSets, universe, hits, minGeneSetSize = 1,
                             pAdjustMethod = "BH", 
                             fix_gene_sets = T,
                             fix_hits = T,
                             fix_universe = T) {
  
  ## Performs a hypergeometric test using the hyperGeoTest function for each gene set in
  ## a collection of gene sets. Automatically sorts output so most significant genes are at the top
  ## of the output.
  ## Inputs:
  ## collectionOfGeneSets - a list of character vectors corresponding to genes in a gene set
  ## universe - a character vector of genes we wish to consider. Genes in the gene sets and hits
  ##            will be subset to these genes
  ## hits - a character vector of genes we consider to be hits
  ## minGeneSetSize - the minimum number of genes in a gene set required for the gene set to be tested
  ## pAdjustmethod - a character vector to be fed to the method argument in p.adjust() for multiple
  ##                 hypothesis correction
  ## Outputs:
  ## matrix with rownames corresponding to the names of the gene sets and columns
  ## as follows:
  ## `Universe Size` - the number of unique genes in the universe
  ## `Gene Set Size` - the number of genes in the tested gene set
  ## `Total hits` - the total number of genes in 'hits'
  ## `Expected Hits` - see hyperGeoTest
  ## `Observed Hits` - see hyperGeoTest
  ## `Odds Ratio` - see hyperGeoTest
  ## `Pvalue` - the pvalue of the hypergeometric test for the corresponding gene set (see hyperGeoTest)
  ## `Adjusted.Pvalue` - the pvalue of the hypergeometric test for the corresponding gene set (Pvalue with
  ##                    multiple hypothesis correction)
  
  if(fix_gene_sets) {
    # Fix gene sets so their elements are unique and in the universe
    gene_sets = lapply(collectionOfGeneSets, function(gene_set) {
      gene_set = unique(gene_set)
      gene_set = intersect(gene_set, universe)
      return(gene_set)
    })
  } else {
    gene_sets = collectionOfGeneSets
  }
  
  # Subset to the gene sets with at least the minimum number of genes
  gene_sets = gene_sets[sapply(gene_sets, length) >= minGeneSetSize]
  
  # Return NULL if there are no gene sets passing our criteria
  if(length(gene_sets) == 0) {
    print("No genesets pass the cutoff size")
    return(NULL)
  }
  
  if(fix_universe) {
    # Make sure there are no repeated genes in the universe
    universe = unique(universe)
  }
  
  if(fix_hits) {
    
    # Make sure there are no repeated genes in the hits
    hits = unique(hits)
    
    # Make sure all the genes in the hits are part of the universe
    hits = intersect(hits, universe)
    
  }
  
  # For each gene set
  dfs = lapply(gene_sets, function(gene_set) {
    # Perform a hypergeometric test on the gene set
    result = hyperGeoTest(gene_set, universe, hits, fix_universe = FALSE, fix_gene_set = FALSE, fix_hits = FALSE)
    # Convert the results from a list to a dataframe
    result = as.data.frame(result)
  })
  
  # Put all the results together into a single data frame
  df = Reduce(rbind, dfs)
  
  # Name the rows of the data frame based on the gene sets
  rownames(df) = names(gene_sets)
  
  # Get the number of genes in each gene set that are in the universe
  gene_set_sizes = sapply(gene_sets, length)
  
  # Format the data frame and add desired information on the sizes of the hits, gene sets
  # and universe
  colnames(df) = c('Pvalue', 'Observed Hits', 'Expected Hits', 'Odds Ratio')
  df$`Universe Size` = length(universe)
  df$`Gene Set Size` = gene_set_sizes
  df$`Total Hits` = length(hits)
  df = df[c('Universe Size', 'Gene Set Size', 'Total Hits', 'Expected Hits',
            'Observed Hits', 'Odds Ratio', 'Pvalue')]
  
  # Add a column for the adjusted pvalue to the data frame
  df$Adjusted.Pvalue = p.adjust(df$Pvalue, method = pAdjustMethod)
  
  # Sort the data frame in increasing order by Pvalue
  df = df[order(df$Pvalue), ]
  
  # Convert the results to a matrix
  result = as.matrix(df)
  
  return(result)
}

multiHyperGeoTests = function(collectionsOfGeneSets, universe, hits, minGeneSetSize = 1,
                              pAdjustMethod = "BH") {
  ## Wrapper function of multiHyperGeoTest for multiple lists of gene sets with convenient formatting
  ## for outputs
  ## Inputs:
  ## collectionsOfGeneSets - list of lists of gene sets, which will be passed to multiHyperGeoTest
  ## universe - same as for multiHyperGeoTest
  ## hits - same as for multiHyperGeoTest
  ## minGeneSetSize - same as for multiHyperGeoTest
  ## pAdjustMethod - same as for multiHyperGeoTest
  ## Outputs:
  ## data frame - analogous to the output of multiHyperGeoTest, but with extra columns:
  ## `source` - the name of the collection from which the gene set came
  ## `across.Adjusted.Pvalue` - multiple hypothesis corrected p value for all gene sets (across collections),
  ##                            whereas Adjusted.Pvalue simply corrects within a collections.
  ## Additionally, this data frame is organized by (within collection) adjust pvalue rather than uncorrected
  ## pvalue as in multiHyperGeoTest
  
  # For each collection of gene sets
  results = mapply(function(collectionOfGeneSets, source) {
    # Run hypergeometric tests on each gene set in the collections
    result = multiHyperGeoTest(collectionOfGeneSets, universe, hits, minGeneSetSize, pAdjustMethod)
    # Convert to a data frame
    result = as.data.frame(result)
    # Add a column for the gene set name
    if(nrow(result) == 0){
      print(paste0("No enrichments could be computed for ", source))
    }else{
      result$source = source
    }
    return(result)
  }, collectionsOfGeneSets, names(collectionsOfGeneSets), SIMPLIFY = FALSE)
  
  # Combine the gene sets from each collection
  results = Reduce(rbind, results)
  # Adjust the pvalues across collections of gene sets
  results$across.Adjusted.Pvalue = p.adjust(results$Pvalue, pAdjustMethod)
  # Sort gene sets in ascending order by Adjusted Pvalue
  if(nrow(results) == 0){
    print("No enrichments could be computed for any genesets")
  }else{
    results = results[order(results$Adjusted.Pvalue, decreasing = FALSE),]
  }
  
  return(results)
}
