get_intersecting_data <- function(eset.list, 
                                  id.col = 'patient_id',
                                  meta.cols = c('condition', 'race', 'gender', "Age")
                                  ){
  require(Biobase)
  require(BiocGenerics)
  # given a list of expressionsets, looks for samples shared between the two and returns
  #  a list of expression matrices and a single metadata dataframe with those samples
  
  # Arguments: eset.list, a list of expression sets
  #        The expression matrices will have the same samples/patients and be in the same order.
  #        id.col is the column that you would like to use to look for intersections
  #        meta cols, character vector, are the metadata columns that you would 
  #          like to appear in the final phenodata matrix

  # value:
  #        a  list with of the expression matrices and one pData matrix
  
  #Get the intersecting id's (patients or samples)
  id.list <- lapply(eset.list, function(eset) eset[[id.col]])
  sample.shared <- Reduce(intersect, id.list)
  
  #Get shared samples and put esets in same order
  eset.list <- lapply(eset.list, function(eset){
    eset <- eset[, match(sample.shared, eset[[id.col]])]
  })
  
  # Make sure the desired id (sample or patient) is used as the column names of the data matrices
  expr.list <- lapply(eset.list, function(eset){
    expr <- exprs(eset)
    colnames(expr) <- eset[[id.col]]
    expr
  })
  
  #Make sure column names are the same between the data matrices
  colnm.list <- lapply(expr.list, colnames)
  colnm.test <- all(sapply(colnm.list[-1], FUN = identical, colnm.list[[1]]))
  stopifnot(colnm.test)
  
  # For the pdat, select the desired columns to keep.
  pdat.list <- lapply(eset.list, function(eset){
    pdat <- pData(eset)
    pdat <- pdat[, c(id.col, meta.cols)]
  })
  # Currently have two pdat's. Make sure that they are in same order
  pdat.test <- all(sapply(pdat.list[-1][[id.col]], FUN = identical, pdat.list[[1]][[id.col]]))
  stopifnot(pdat.test)
  
  #Select first pdat as the output pdat
  pdat.shared <- pdat.list[[1]]
  
  out <- list(expr = expr.list, pdat = pdat.shared)
  
  out
}
