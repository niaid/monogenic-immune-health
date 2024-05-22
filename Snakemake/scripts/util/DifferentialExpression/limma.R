#https://support.bioconductor.org/p/21333/
#https://kasperdanielhansen.github.io/genbioconductor/html/limma.html
#http://jtleek.com/genstats/inst/doc/03_09_Calculating_statistics.html#fit-many-statistics-with-limma
#http://lgatto.github.io/tipsntRicks/bioconductor/2012/05/13/about-limma/
#https://support.bioconductor.org/p/67395/

library(dplyr)

# Source utility file to get the AI, PID, and TERT-TERC condition groups
source('scripts/util/Groups/groups.R')

# Set a global containing the set of all conditions in the training data
CONDITIONS = util.get_conditions()

# Set a global containing the set of all AI conditions
AI_CONDITIONS = util.get_ai()

# Set a global containing the set of all PID conditions
PID_CONDITIONS = util.get_pid()

# Set a global containing the set of all TERT-TERC conditions
TERT_TERC_CONDITIONS = util.get_tert_terc() 

add_meds_eset <- function(eset, med.dat){
  #adds medication data to the expressionset
  
  #remove the patient_id column from med.dat
  med.dat = med.dat[,setdiff(colnames(med.dat), 'patient_id')]
  
  #add medication columns to expression set
  new.pheno <- left_join(pData(eset), med.dat, by = "visit_id")
  phenoData(eset) <- AnnotatedDataFrame(new.pheno) 
  
  eset
}

make_design <- function(eset, group.col = "condition", ref.level = "Healthy"){
  ## Makes a design matrix from a sample-level eset using a reference level of healthy patients, 
  ## with information on the patient ids, condition groups, sex, and age.
  ## Inputs:
  ## eset - a sample level eset
  ## group.col - the name of the column in pData(eset) with the condition labels
  ## ref.level - the reference condition for which the model should be created
  ## Output:
  ## limma design matrix - a design matrix that can be fit using limma
  
  # Get the condition group of each sample
  group <- factor(eset[[group.col]])
  
  # If the reference level is provided, make it the factor reference level
  if(!is.null(ref.level)){
    group <- relevel(group, ref.level)
  }
  
  # Get the sexes and ages associated with each sample
  sex <- factor(eset$gender)
  age <- eset$patient_age_at_time_of_blood_draw
  visit_type <- eset$visit_type
  
  design <- model.matrix(~group + age + sex + visit_type) #cannot include batch as it makes matrix linearly dependent as error in duplicateCorrelation step below (https://support.bioconductor.org/p/32455/)
  
  # Reformat the intercept column
  colnames(design)[colnames(design) == "(Intercept)"] <- "Intercept"
  
  # Ensure that the column names of the design matrix are all proper names
  colnames(design) <- make.names(colnames(design))
  
  design
}

make_sex_linked_design <- function(eset, 
                        group.col = "condition",
                        sex.col = "gender",
                        group.ref.level = "Healthy",
                        sex.ref.level = "M") {
  ## Makes a design matrix from a sample-level eset using a reference level of healthy patients, 
  ## with information on the patient ids, condition groups, sex, and age. Includes a term for the interaction
  ## between sex and condition
  ## Inputs:
  ## eset - a sample level eset
  ## group.col - the name of the column in pData(eset) with the condition labels
  ## group.ref.level - the reference condition for which the model should be created
  ## sex.ref.level - the reference sex for which the model should be created
  ## Output:
  ## limma design matrix - a design matrix that can be fit using limma
  
  # Get the condition group of each sample
  group = factor(eset[[group.col]])
  sex = factor(eset[[sex.col]])
  
  # If the reference level is provided, make it the factor reference level
  if(!is.null(group.ref.level)) {
    group = relevel(group, group.ref.level)
  }
  if(!is.null(sex.ref.level)) {
    sex = relevel(sex, sex.ref.level)
  }
  
  # Get the sexes and ages associated with each sample
  age = eset$patient_age_at_time_of_blood_draw
  visit_type = eset$visit_type
  
  design <- model.matrix(~group*sex + age + visit_type)
  
  # Reformat the intercept column
  colnames(design)[colnames(design) == "(Intercept)"] <- "Intercept"
  
  # Ensure that the column names of the design matrix are all proper names
  colnames(design) <- make.names(colnames(design))
  
  design
}

fit_limma <- function(eset, design){
  ## Fit a design matrix using data from a sample level expression set
  ## Inputs:
  ## eset - the sample-level eset whose pData was used to create the design matrix
  ## design - the design matrix made using the make_design function
  ## Outputs:
  ## fit - a limma fit object from fitting the design matrix
  
  # Create a correlation object using duplicate correlations to turn the patient_id
  # into a random effect
  corfit <- duplicateCorrelation(eset, design, block = factor(eset$patient_id))
  # Fit the model to the data in the eset using the design and blocking by patient
  fit <- lmFit(eset, design, block = eset$patient_id, correlation = corfit$consensus)
  fit
}

make_condition_groups = function(conditions) {
  ## Make a list containing a group for each condition in the conditions vector
  ## that can be used to easily create contrast matrices using the make_contrasts_mat function
  ## Inputs:
  ## conditions - character vector of all conditions to be included
  ## Outputs:
  ## list - condition groups (with one condition for each group).
  ## Groups and conditions have the same name.
  
  # Ensure there are no duplicate conditions
  conditions = unique(conditions)
  conditions = sort(conditions)
  # Turn each condition into a group
  conditions = as.list(conditions)
  # Name the groups the same as the conditions
  names(conditions) = conditions
  return(conditions)
}

make_AI_PID_groups = function(ai.conditions, 
                              pid.conditions,
                              tert.terc.conditions) {
  ## Make a list containing groups of AI, PID, and TERT-TERC conditions
  ## Inputs:
  ## ai.conditions - a character vector with the AI conditions
  ## pid.conditions - a character vector with the PID conditions
  ## tert.terc.conditions - a character vector with the TERT-TERC conditions
  ## Outputs:
  ## list - list of AI, PID, and TERT-TERC condition groups
  
  # get all AI patients
  ai.conditions = unique(ai.conditions)
  ai.conditions = sort(ai.conditions)
  
  # Get all PID conditions
  pid.conditions = unique(pid.conditions)
  pid.conditions = sort(pid.conditions)
  
  # Get all TERT-TERC conditions 
  tert.terc.conditions = unique(tert.terc.conditions)
  tert.terc.conditions = sort(tert.terc.conditions)
  
  # Create a list of the condition groups
  list(AI = ai.conditions, PID = pid.conditions, TERT.TERC = tert.terc.conditions)
}

make_contrasts_mat = function(fit, groups, cross = FALSE) {
  ## Make a constrasts matrix with variables corresponding to the fit. The
  ## groups decribe the condition groups to be used in the contrast matrix.
  ## If 'cross' is FALSE, the contrast matrix is made in relation to healthy. 
  ## If 'cross' is TRUE, the contrast matrix is made in relation to all other conditions 
  ## across any of the groups in 'groups'.
  ## Inputs:
  ## fit - A fit object as derived using fit_limma
  ## groups - a named list of condition groups
  ## cross - FALSE for versus healthy; TRUE for versus all
  
  # Get the covariates used in the fit
  vars = colnames(fit)
  
  # Subset the groups to just those from the fit
  groups = lapply(groups, function(group) {
    select = make.names(paste0('group', group)) %in% vars
    group = group[select]
  })
  
  # Remove any empty groups
  groups = groups[sapply(groups, length) > 0]
  
  # Get the set of all conditions in the groups
  all.conditions = unname(unique(unlist(groups, recursive = T)))
  
  # For each group
  mat = sapply(groups, function(group) {
    
    # Find the set of conditions in groups not in this group
    background = setdiff(all.conditions, group)
    
    # Create an empty vector with one entry per covariate
    x = rep(0, length(vars))
    names(x) = vars
    
    # Convert the condition names to the same format as the covariates
    group = make.names(paste0('group', group))
    # Ensure that all the conditions in the group were included as covariates in the fit
    stopifnot(all(group %in% vars))
    # Split the positive contrast weight between the conditions in the group
    x[group] = 1 / length(group)
    
    # If we wish to put the contrasts matrix in the versus all format
    if(cross) {
      # Ensure the background group names are of the correct format
      background = make.names(paste0('group', background))
      # Ensure that all the conditions in the background were included as covariates in the fit
      stopifnot(all(background %in% vars))
      # Divide the negative contrast weight between the background conditions
      x[background] = -1 / length(background)
    }
    
    return(x)
  })
  
  colnames(mat) = names(groups)
  return(mat)
}

get_topTable_list <- function(fit) {
  #Performs T test for every group and returns list of dataframe with associated info
  
  # test.names refers to the contrast/group that is being tested
  test.names = colnames(fit$coefficients)
  
  #Ierate through each contrast/group and use limma::topTable to compute t-statistics and p values
  toptab.list <- lapply(test.names, function(test.name){
    x=topTable(fit, coef=test.name, adjust.method = "BH", sort.by = "none", number = nrow(fit[["p.value"]]))
  })
  
  names(toptab.list) = test.names
  
  toptab.list
}

toptab_list_to_matlist <- function(toptab.list){
  # converts list of dataframes for each group to feature X group matrices
  # with t stats and p values
  
  #Make sure that all of the rownames are the same for all topTables/contrasts
  # Set one toptable as reference and make sure all match
  feat <- rownames(toptab.list[[1]]) #setting first toptable rownames as reference 
  test <- sapply(toptab.list, function(x) {
    identical(feat, rownames(x))
  })
  stopifnot(all(test))
  
  #Make t statistic matrix for all contrasts
  t <- sapply(toptab.list, function(x) x[["t"]])
  #Make p value matrix for all contrasts
  P.Value <- sapply(toptab.list, function(x) x[["P.Value"]])
  #Make fdr adjust p value matrix for all contrasts
  adj.P.Val <- sapply(toptab.list, function(x) x[["adj.P.Val"]])
  
  #Put matrices into list as output
  out.list <- list(t = t, P.Value = P.Value, adj.P.Val = adj.P.Val)
  
  #Make sure all matrices have the featurenames as rownames
  out.list <- lapply(out.list, function(x){
    rownames(x) <- feat
    x
  })
  
  out.list
}

get_ebayes_stats = function(fit) {
  # Wrapper function for the ebayes steps
  fit.ebayes = eBayes(fit)
  toptab.list = get_topTable_list(fit.ebayes)
  ebayes.t.mat.list = toptab_list_to_matlist(toptab.list)
  ebayes.t.mat.list$effect.size = fit$coefficients
  ebayes.t.mat.list$cohens.d = fit$coefficients / sqrt(fit.ebayes$s2.post)
  
  ebayes.t.mat.list
}

get_traditional_stats = function(fit) {
  # returns feature X group matrices
  # performs a standard t-test without empirical bayes moderation
  # see https://support.bioconductor.org/p/113833/
  
  # Compute t statistics and p values- objects are feature X contrast matrices
  t.stat <- apply(fit$coefficients/fit$stdev.unscaled, 2, function(x){x/fit$sigma})
  t.stat.p.value <- 2 * pt(-abs(t.stat), df = fit$df.residual)
  t.stat.adj.p <- apply(t.stat.p.value, 2, p.adjust, method = "BH")
  
  #outputs list of matrices
  out = list(
    t = t.stat,
    P.Value = t.stat.p.value,
    adj.P.Val = t.stat.adj.p,
    effect.size = fit$coefficients,
    cohens.d = fit$coefficients / sqrt(fit$sigma)
  )
  
  out
}
