util.make.signature = function(score, mat, n = 150, cor = .35, method = 's') {
  ## Makes a surrogate for the vector 'score' using features from the matrix 'mat'.
  ## The signature is composed of the features in 'mat' most correlated with 'score', with a cap of 'n'
  ## features, and a minimum correlation threshold (calculated using method 'method') of 'cor'.
  ## The signature is divided into the features that are positively correlated with 'score'
  ## and the features that are negatively correlated with 'score'
  ## Inputs:
  ## score - a numeric vector of scores for each subject
  ## mat - a numeric matrix of dimensions subjects x features, where the rows
  ##       are in the same order as the entries of 'score'
  ## n - the maximum number of features in the signature
  ## cor - the minimum (absolute) value of the correlation to be included in the
  ##       signature
  ## method - a string specifying the method used to calculate the correlations
  
  # Get all gene correlations associated with a score
  feature.cors = cor(score, mat, method = method)
  
  # Convert the correlations from a matrix to a vector
  x = feature.cors[1,]
  names(x) = colnames(feature.cors)

  # Rearrange the vector so the correlations are in decescending (absolute order)
  x = x[order(abs(x), decreasing = TRUE)]
  
  # Extract the features from the top n that have a correlation greater than the threshold
  positive = intersect(names(x)[x >= cor], names(x)[1:n])
  
  # Extract the features from the top n that have a correlation less than the negative threshold
  negative = intersect(names(x)[x <= -cor], names(x)[1:n])
  
  # Make the signature
  signature = list(positive = positive, negative = negative)
  
  return(signature)
}