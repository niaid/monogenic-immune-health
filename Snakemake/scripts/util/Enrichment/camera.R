library(limma)

make_indices <- function(geneset.list, universe, min.geneset.size){
  # Get indices object that can be used for cameraPR
  
  # Arguments:
  #  geneset.list: a list containing character vectors of genesets
  #  universe: the background, all genes that should be used as background for enrichment,
  #    typically the rownames of the expressionset
  #  min.geneset.size: The smallest number of genes allowed in a geneset,
  #    genesets with fewer genes than min.geneset.size in the universe will be removed
  
  # Value
  #   Returns a list of vectors of indices that can be used by limma::cameraPR to compute 
  
  indices <- lapply(geneset.list, function(geneset){
    which(universe %in% geneset)
  })
  indices <- indices[sapply(indices, length) >= min.geneset.size]
  
  indices
}

get_cormat <- function(x, y, meth = "spearman"){
  # reorder matrices so that samples match and compute correlation matrix
  # output can be used in cameraPR_cor
  
  # arguments:
  #   x, y: matrices, samples X features
  #   meth: the method to compute the correlation
  
  # value
  #   numeric matrix with correlation values
  
  intersection <- intersect(rownames(x), rownames(y))
  x <- x[match(intersection, rownames(x)), ]
  y <- y[match(intersection, rownames(y)), ]
  
  cormat <- cor(x, y, method = meth)
  cormat
}

cameraPR_cor <- function(cormat, geneset.list, min.geneset.size, use.ranks = FALSE, abs.cor = FALSE){
  #Takes a correlation matrix as input and performs limma's cameraPR gene set test
  
  #Arguments:
  #  cormat: correlation matrix that can be computed with get_cormat
  #    Features being tested for enrichment are in rows. Genes are in columns of the correlation matrix
  #    rows and columns of the correlation matrix must have names
  #  geneset.list: a list of character vectors with the gene names for each set
  #  min.geneset.size: The smallest number of genes allowed in a geneset,
  #    genesets with fewer genes than min.geneset.size in the universe will be removed
  #  use.ranks: 
  #    TRUE: perform cameraPR with a wilcoxon
  #    FALSE: perform cameraPR with a t-test
  #  abs.cor: boolean, take the absolute value of the correlation?
  #
  #Value:
  #  a list of dataframes output from cameraPR 
  #  each dataframe corresponds to a single row of the input correlation matrix with the same name
  
  universe <- colnames(cormat)
  indices <- make_indices(geneset.list, universe, min.geneset.size)
  
  if(abs.cor){
    cormat <- abs(cormat)
  }
  
  enrich.list <- lapply(rownames(cormat), function(row.nm){
    corvec <- cormat[row.nm, , drop = TRUE]
    
    dat <- cameraPR(corvec, indices, use.ranks, sort = FALSE, inter.gene.cor = 0)
    dat$geneset <- names(indices)
    dat
  })
  names(enrich.list) <- rownames(cormat)
  
  enrich.list
}

multiGST <- function(indices, statistics, alternative, type = "auto"){
  #Not currently used in monogenic pipeline
  
  # see limma gene set test documentation to understand these arguments.
  # This is simply a wrapper that outputs multiple test results as a dataframe
  
  p.vals <- vector("numeric", length = length(indices))
  names(p.vals) <- names(indices)
  
  for(gset in names(p.vals)){
    index <- indices[[gset]]
    p.vals[[gset]] <- geneSetTest(index, statistics, alternative = alternative, 
                                  type= type, ranks.only = TRUE)
  }
  
  out.dat <- data.frame(geneset = names(p.vals),
                        NGenes = sapply(indices, length),
                        PValue = p.vals,
                        FDR = p.adjust(p.vals, method = "fdr"))
  
  out.dat
  
}

multiGST_cor <- function(cormat, geneset.list, min.geneset.size, type = "auto", alternative, use.ranks = TRUE, abs.cor, cor.meth){
  #Not currently used in monogenic pipeline
  
  #Takes a correlation matrix as input and performs limma's cameraPR gene set test
  
  universe <- colnames(cormat)
  indices <- make_indices(geneset.list, universe, min.geneset.size)
  
  if(abs.cor){
    cormat <- abs(cormat)
  }
  
  enrich.list <- lapply(rownames(cormat), function(row.nm){
    rho <- cormat[row.nm, , drop = TRUE]
    
    dat <- multiGST(indices, rho, alternative, type = type)
    dat
  })
  names(enrich.list) <- rownames(cormat)
  
  enrich.list
}


