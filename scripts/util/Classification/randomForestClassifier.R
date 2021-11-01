library(randomForest)

cross.validation = function(X, y, pos = "Healthy") {
  ## Runs LOO cross validations using a random forest classifier
  ## on a data set and returns the LOO CV probability of being the sample
  ## being from the positive class (specified by the argument 'pos')
  ## Inputs:
  ## X - a numeric design matrix of dimensions samples x features
  ## y - the response vector
  ## pos - the set of values of the response vector considered to be
  ##       the positive class
  ## Outputs:
  ## The LOO CV probabilities of being part of the positive class for each sample
  
  # Convert the response vector values to 0s and 1s (negative and positive classes)
  y = ifelse(y %in% pos, 1, 0)
  # Change y to a factor
  y = factor(y, levels = c(0,1))
  
  # For each index
  predictions = sapply(1:length(y), function(test.index) {
    # Create a training matrix without that index
    X.train = X[-test.index, , drop = FALSE]
    # Get the training samples' response values
    y.train = y[-test.index]
    # Create a testing matrix with only that index
    X.test = X[test.index, , drop = FALSE]
    # Create a model using the training data
    model = randomForest(X.train, y.train)
    # Predict the value of the testing set
    probabilities = predict(model, X.test, type = 'prob')
    # Return the probability that the sample is of the positive class
    probabilities[,'1']
  })
  
  return(predictions)
}

get.rf.model = function(X, y, pos = 'Healthy') {
  ## Trains a random forest classifier and return the trained model
  ## Inputs:
  ## X - a numeric design matrix of dimensions samples x features
  ## y - the response vector
  ## pos - the set of values of the response vector considered to be
  ##       the positive class
  ## Outputs:
  ## The random forest model trained on the inputs and outputs
  
  # Convert the response vector values to 0s and 1s (negative and positive classes)
  y = ifelse(y %in% pos, 1, 0)
  # Change y to a factor
  y = factor(y, levels = c(0,1))
  
  # Train a random forest model
  model = randomForest(X, y)
  
  return(model)
}

get.gvis = function(X, y, pos = 'Healthy') {
  ## Trains a random forest classifier and extracts the GVIs for each features
  ## Inputs:
  ## X - a numeric design matrix of dimensions samples x features
  ## y - the response vector
  ## pos - the set of values of the response vector considered to be
  ##       the positive class
  ## Outputs:
  ## The GVIs associated with each feature in the derived random forest

  # Train a random forest model
  model = get.rf.model(X, y, pos = pos)
  
  # Extract the gvis from the trained model
  model$importance[,'MeanDecreaseGini']
}