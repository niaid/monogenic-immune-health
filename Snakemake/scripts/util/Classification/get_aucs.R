require(ggplot2)

get_roc = function(x, y, pos = 1) {
  ## Constructs an ROC curve based on score x with binary response vector y
  ## Inputs:
  ## x - numeric vector corresponding to likelihoods of associated samples being of the positive class
  ## y - binary response vector or factor corresponding to the samples in x
  ## pos - the set of possible values for y that should be associated with the positive class
  ## Outputs:
  ## data frame with false positive ratio and true positive ratio pairs on the ROC curve
  
  # Convert y to a vector of TRUEs and FALSEs corresponding to if y is in the positive class
  y = y %in% pos
  
  # Sort the scores so they are in descending order
  order = order(x, decreasing = TRUE)
  x = x[order]
  # Change the response vector so the samples correspond with those in x
  y = y[order]
  
  # Get the total number of samples
  n = length(y)
  # Get the total number of positive samples
  P = sum(y)
  # Get the total number of negative samples
  N = n - P
  
  # Initialize vectors of the true positive ratios and false positive ratios
  tprs = 0
  fprs = 0
  for(i in 1:n) {
    y_subset = y[1:i]
    # Get the number of true positives in the first i samples
    TP = sum(y_subset)
    # Get the number of false positives in the first i samples
    FP = i - TP
    # Get the true positive ratio if we were to guess the first i samples
    # were positive
    tpr = TP / P
    # Get the false positive ratio if we were to guess the first i samples
    # were positive
    fpr = FP / N
    # Append the true positive ratio to our list
    tprs = c(tprs, tpr)
    # Append the false positive ratio to our list
    fprs = c(fprs, fpr)
  }
  # Return a data frame with the tprs and fprs
  roc = data.frame(tpr = tprs, fpr = fprs)
  
  return(roc)
}

get_auc = function(roc) {
  ## Gets the auc from an roc curve
  ## Inputs:
  ## roc - Data frame corresponding to an roc curve with columns $fpr and $tpr,
  ##       corresponding to false positive ratios and true positive ratios along
  ##       the auc curve
  ## Outputs:
  ## auc - the area under the ROC curve
  
  # Initialize the auc to 0
  auc = 0
  # For each step along the curve 
  for(i in 2:nrow(roc)) {
    # Add the area associated with that step to the auc
    auc = auc + (roc$fpr[i] - roc$fpr[i-1]) * roc$tpr[i]
  }
  return(auc)
}

get_auc_pvalue = function(x, y, pos = 1) {
  x = seq(0, 1, .1)
  y = c(0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1)
  roc = get_roc(x, y, pos = 1)
  auc = get_auc(roc = roc)
  
  auc.perms = sapply(1:1000, function(i) {
    y.perm = sample(y)
    roc = get_roc(x, y.perm)
    auc = get_auc(roc)
  })
  
  p = mean(auc < auc.perms)
}