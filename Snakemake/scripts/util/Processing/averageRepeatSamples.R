library(Biobase)

averageRepeatSamples = function(samples.eset,
                                patient.id.col = 'patient_id',
                                meta.cols = c('condition', 'race', 'gender'), 
                                age.col = 'patient_age_at_time_of_blood_draw') {
  
  ## Utility function to average data (and ages) over samples from the same patient
  ## Inputs:
  ## samples.eset - an eset with data at the sample level
  ## patient.id.col - the name of the column in fData(samples.eset) that identifies patient ids
  ## meta.cols - the other names of the columns in fData(samples.eset) to be
  ##             included in the feature data of the sample level eset. Can optionally be empty
  ## age.col - the name of the age column in fData(samples.eset) to be averaged over across samples from
  ##           the same patient. If NA, the column will be left out
  ## Outputs:
  ## an eset with data from samples.eset, but averaged to be at the subject level
  
  # Extract the data matrix from the eset
  samples.matrix = exprs(samples.eset)
  
  # Add the ages as a feature (if age.col is not NA)
  if(!is.na(age.col)) {
    # Get the ages associated with each sample
    samples.ages = pData(samples.eset)[[age.col]]
    # Convert the sample ages from a matrix to a vector
    samples.ages = matrix(samples.ages, ncol = ncol(samples.matrix))
    # Name the sample ages with the visit names
    colnames(samples.ages) = colnames(samples.matrix)
    # Name the sample ages row
    rownames(samples.ages) = 'Age'
    # Combine the sample ages with the data matrix
    samples.matrix = rbind(samples.matrix, samples.ages)
  }
  
  # Split the columns into groups based on the patient ids
  sample.groups = split(colnames(samples.matrix), samples.eset[[patient.id.col]])
  
  # Get the mean of each group of columns corresponding to repeat samples and put them into a matrix
  subjects.matrix = sapply(sample.groups, function(x) { rowMeans(samples.matrix[, x, drop = FALSE]) })
  # Name the columns with the patient ids
  colnames(subjects.matrix) = names(sample.groups)
  
  # Extract the sample and feature meta data from the eset
  samples.meta = pData(samples.eset)
  features.meta = fData(samples.eset)
  
  ## Here we ensure that desired metadata corresponding to each 
  ## group of samples is identical
  
  # Subset sample level meta data to desired columns
  check.meta = samples.meta[, c(patient.id.col, meta.cols), drop = FALSE]
  # Get groups of samples for checking
  check.groups = split(rownames(check.meta), check.meta[[patient.id.col]])
  # For each group of samples
  check.meta = sapply(check.groups, function(x) {
    # For each meta data feature, get the number of unique values 
    apply(check.meta[x, , drop = FALSE], 2, function(y) {length(unique(y))}) 
  })
  # If there is more than one value for any column, throw an error
  stopifnot(all(check.meta == 1))
  
  ## Done with check ##
  
  # Convert the meta data at the sample level into meta data at the subject level
  subjects.meta = samples.meta[! duplicated(samples.meta[[patient.id.col]]), c(patient.id.col, meta.cols), drop = FALSE]
  # Name the subject level meta data with the patient ids
  rownames(subjects.meta) = subjects.meta[[patient.id.col]]
  # Ensure that the subject level data and subject level meta data are in the same order
  subjects.meta = subjects.meta[names(sample.groups), , drop = FALSE]
  
  # Extract and remove the ages as a feature (if age.col is not NA)
  if(!is.na(age.col)) {
    # Get the mean age of each patient
    subjects.ages = subjects.matrix['Age',]
    # Add these ages as the ages column in the metadata
    subjects.meta$Age = subjects.ages
    # Remove the age from the data matrix
    subjects.matrix = subjects.matrix[setdiff(rownames(subjects.matrix), 'Age'), , drop = FALSE]
  }
  
  # Ensure that the data and feature data are in the same order
  features.meta = features.meta[rownames(subjects.matrix), , drop = FALSE]
  
  # Assemble the expression set
  subjects.eset = ExpressionSet(subjects.matrix)
  phenoData(subjects.eset) = AnnotatedDataFrame(subjects.meta)
  featureData(subjects.eset) = AnnotatedDataFrame(features.meta)
  
  return(subjects.eset)
}