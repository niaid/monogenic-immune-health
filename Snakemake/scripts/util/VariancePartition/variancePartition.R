library(variancePartition)

patient_id_variance_partition = function(eset, patient.id.col = 'patient_id') {
  ## Run a variance partition where the only covariate is the patient
  ## Inputs:
  ## eset - an expression set with the parameters to study using
  ##        the variance partition, and metdata containing the patient_id
  ## patient.id.col - the name of the column in pData(eset) containing the patient_ids
  ## Output:
  ## variance partition object containing estimated variance explained by the patient.id
  ## random effect and the residual variances for each parameter
  
  # Make the sample information matrix
  info = pData(eset)[patient.id.col]
  colnames(info) = c('Patient')
  info$Patient = factor(paste('Patient', '_', info$Patient))
  
  # List the fixed and random effects (only covariate is patient)
  form = ~ (1|Patient)
  
  # Fit the model
  varPart <- fitExtractVarPartModel(exprs(eset), form, info)
  
  # Clean the results
  vp <- sortCols( varPart )
  
  return(vp)
  
}

condition_medication_variance_partition = function(eset, medications) {
  ## Runs variance partition across samples with covariates for patients, condition, and several
  ## medication groups. All medications passed in the 'medications' matrix are used.
  ## Inputs:
  ## eset - the sample-level Expression Set on which to perform the variance partition.
  ## medications - a binary matrix outlining which medications a patient is on, with rows in the medication
  ##               matrix corresponding to columns of the eset 
  ## Outputs:
  ## a variance partition object estimating the amount of variance that can be assigned to each
  ## covariate
  
  # Get the sample meta inforamtion (patient ids and conditions)
  info = pData(eset)[c('patient_id', 'condition')]
  # Rename the columns of info
  colnames(info) = c('Patient', 'Condition')
  # Initialize the formula for the variance partiton
  form = '~ (1|Patient) + (1|Condition)'
  
  # Combine the medication information with the patient ids and conditions
  info = cbind(info, as.data.frame(medications))
  # Add each medication as a random effect in the formula
  for(medication in colnames(medications)) {
    form = paste0(form, ' + (1|', medication,')')
  }
  # Print the formula so the user can confirm the covariates they wish for are present 
  print(form)
  
  # Run the variance partition model
  varPart = fitExtractVarPartModel(exprs(eset), form, info)
  vp = sortCols( varPart )
  
}