# Load the monogenic database
load('Inputs/Metadata/monogenic.de-identified.metadata.RData')


util.get_ai = function() {
  ## Gets AI groups
  ## Outputs:
  ## character vector - all the conditions belonging to the AI group
  monogenic.all.assays = monogenic.all.assays[monogenic.all.assays$analysis_group == 'Discovery', ]
  conditions = monogenic.all.assays$condition[monogenic.all.assays$condition_w_autoinflammatory == "yes"]
  conditions = unique(conditions[!is.na(conditions)])
  return(conditions)
}

util.get_tert_terc = function() {
  ## Gets TERT-TERC groups
  ## Outputs:
  ## character vector - all the conditions belonging to the TERT-TERC group
  return(c('TERT', 'TERC'))
}

util.get_pid = function() {
  ## Gets PID groups. Note that we define the PID group not to include DADA2 or the NEMO carrier
  ## Outputs:
  ## character vector - all the conditions belonging to the PID group
  monogenic.all.assays = monogenic.all.assays[monogenic.all.assays$analysis_group == 'Discovery', ]
  conditions = monogenic.all.assays$condition[monogenic.all.assays$condition_w_immunodeficiency == "yes"]
  conditions = unique(conditions[!is.na(conditions)])
  others = c(util.get_ai(), util.get_tert_terc(), 'NEMO carrier')
  conditions = setdiff(conditions, others)
  return(conditions)
}

util.get_conditions = function() {
  ## Get all conditions seen in the training samples
  ## Outputs:
  ## character vector - all the conditions belonging to training samples (not including Healthy)
  monogenic.all.assays = monogenic.all.assays[monogenic.all.assays$analysis_group == 'Discovery', ]
  conditions = unique(monogenic.all.assays$condition)
  conditions = setdiff(conditions, 'Healthy')
  return(conditions)
}
