library(Biobase)

averageTechnicalReplicates = function(replicates.eset,
                                      visit.id.col = 'visit_id',
                                      meta.cols = c('patient_id',
                                                    'condition',
                                                    'race',
                                                    'gender',
                                                    'assay_desc',
                                                    'patient_age_at_time_of_blood_draw')) {
  
  ## Average data over technical replicates and save selected metadata columns
  ## Note: here we've assumed patient ages, if included, are the same for all samples with the same visit_ids,
  ## and thus do not need to be averaged over
  
  # Extract the replicate level data matrix from the eset
  replicates.matrix = exprs(replicates.eset)
  
  # Split the replicate level data matrix columns into groups based on the corresponding visit ids
  replicate.groups = split(colnames(replicates.matrix), replicates.eset[[visit.id.col]])
  
  # Get the mean of each group of replicates corresponding to the same visit id and put them into a sample level matrix
  samples.matrix = sapply(replicate.groups, function(x) { rowMeans(replicates.matrix[, x, drop = FALSE]) })
  # Name the columns of the sample level matrix as the corresponding visit ids 
  colnames(samples.matrix) = names(replicate.groups)
  
  # Extract the phenotype and feature meta data from the replicate level eset
  replicates.meta = pData(replicates.eset)
  features.meta = fData(replicates.eset)
  
  ## Here we ensure that desired metadata corresponding to each 
  ## group of technical replicates is identical
  
  # Subset replicate level meta data to desired columns
  check.meta = replicates.meta[, c(visit.id.col, meta.cols), drop = FALSE]
  # Get groups of replicates for checking
  check.groups = split(rownames(check.meta), check.meta[[visit.id.col]])
  # For each group of replicates
  check.meta = sapply(check.groups, function(x) {
    # For each meta data feature, get the number of unique values 
    apply(check.meta[x, , drop = FALSE], 2, function(y) {length(unique(y))}) 
  })
  # If there is more than one value for any column, throw an error
  stopifnot(all(check.meta == 1))
  
  ## Done with check ##
  
  # Convert the meta data at the replicate level into meta data at the sample level
  samples.meta = replicates.meta[! duplicated(replicates.meta[[visit.id.col]]), c(visit.id.col, meta.cols), drop = FALSE]
  # Name the rows of the sample level meta data based on the visit ids
  rownames(samples.meta) = samples.meta[[visit.id.col]]
  # Ensure the sample level phenotype data corresponds to the sample level data matrix
  samples.meta = samples.meta[names(replicate.groups), , drop = FALSE]
  
  # Assemble the expression set
  samples.eset = ExpressionSet(samples.matrix)
  phenoData(samples.eset) = AnnotatedDataFrame(samples.meta)
  featureData(samples.eset) = AnnotatedDataFrame(features.meta)
  
  return(samples.eset)
}