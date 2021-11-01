# Base function: gets average z score over all features in a signature, for each patient
## Inputs: X - a matrix with dimensions: subjects x features
##         s - a named list with the following structure:
##             List of 2
##              $ positive: chr - the names of the genes used as positive correlates in the signature
##              $ negative: chr - the name of the genes used as negative correlates in the signature
## Ouptut: A vector with signature scores for each subject
util.get_signature_score = function(X, s) {
  # Get the Z scores of each feature for each patient
  Z = scale(X)
  
  # Only inlcude features from the signature that are in the data
  s.plus = intersect(s$positive, colnames(Z))
  s.minus = intersect(s$negative, colnames(Z))
  
  # Make a matrix with positive correlates if there are any
  if(length(s.plus) > 0) {
    Z.plus = Z[, s.plus, drop = FALSE]
  } else {
    Z.plus = NULL
  } 
  
  # Make a matrix with the negative correlates if there are any
  if(length(s.minus) > 0) {
    Z.minus = -Z[, s.minus, drop = FALSE]
  } else {
    Z.minus = NULL
  }
  
  # If the signature is empty, return NAs
  if(length(s.plus) + length(s.minus) == 0) {
    signature_score = rep(NA, nrow(Z))
    names(signature_score) = rownames(Z)
    return(signature_score)
  }
  
  # Otherwise return the average z-score (after sign-adjustment)
  Z = cbind(Z.plus, Z.minus)
  signature_score = rowMeans(Z)
  return(signature_score)
}

# Wrapper function: gets average z score over all features in a signature, for each patient and signature
## Inputs: X - a matrix with dimensions subjects x features
##         signatures - a named list of signatures (as described above)
## Ouput: Matrix with the signature scores of each subject (of dimension subject x signatures)
util.get_signature_scores = function(X, signatures) {
  signature_scores = sapply(signatures, function(s) {util.get_signature_score(X, s)})
  colnames(signature_scores) = names(signatures)
  return(signature_scores)
}

