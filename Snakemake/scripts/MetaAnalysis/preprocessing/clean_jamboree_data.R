## Organizes the jamboree data together into a metaintegrator analysis object

# Load libraries
library(MetaIntegrator)
library(preprocessCore)

# Set globals
## Compiled jamboree data
DATA.IN.PATH = snakemake@input[[1]]#'Reference/jamboree/processed/jamboree_data.RDS'
## Clean cgps info
CGPS.IN.PATH = snakemake@input[[2]]#'Reference/jamboree/data_analysis_ready/cgps_clean.RDS'

## Cleaned and compiled meta-integrator object
META.OUT.PATH = snakemake@output[[1]]#'Reference/jamboree/data_analysis_ready/meta_studies.RDS'

# Load data
data.sets = readRDS(DATA.IN.PATH)
cgps = readRDS(CGPS.IN.PATH)

# Remove SLE as we do not plan to use these studies due to the small number of covered genes
data.sets = data.sets[names(data.sets) != 'SLE']

# Change data frames to matrices
data.sets = lapply(data.sets, as.matrix)

# Get the genes each disease's data set has in common
genes = Reduce(intersect, lapply(data.sets, rownames))

# Put together metaObject
## Create an empty meta object
metaObj = list()
## For each disease
for(disease in names(data.sets)) {
  ## Get all the studies in this diseass
  studies = names(cgps[[disease]])
  ## Get the set of data corresponding to this disease
  X = data.sets[[disease]]
  ## For each study
  for(study in studies) {
    ## Get the cases and controls associated with this study
    cgp = cgps[[disease]][[study]]
    case.gsms = cgp$cases
    control.gsms = cgp$controls
    gsms = c(case.gsms, control.gsms)
    ## Create an empty object to hold the data
    dataObj = list()
    
    ## Create a data frame with the gsms for the pheno slot in the dataObj
    dataObj$pheno = data.frame(gsms = gsms)
    rownames(dataObj$pheno) = gsms
    
    ## Set the classes associated with each case to be 1 and 
    ## controls to be 0
    dataObj$class = ifelse(gsms %in% case.gsms, 1, 0)
    ## Name the classes with the gsms
    names(dataObj$class) = gsms
    
    ## Set the name associated with the study in the data object to be the study name
    dataObj$formattedName = study
    
    ## Get the data for the desired genes and gsms in this study
    expr = X[genes, gsms]
    ## Remove any subjects with NAs
    expr = expr[rowSums(is.na(expr))==0,]
    ## Normalize quantiles within the study, as if some studies were already quantile normalized,
    ## all of them should be
    expr.normalized = normalize.quantiles(expr)
    ## Ensure all values are positive
    if(min(expr.normalized) <= 0) {
      expr.normalized = expr.normalized - min(expr.normalized) + 1
    }
    ## Name the normalized expression matrix the same as the original
    colnames(expr.normalized) = colnames(expr)
    rownames(expr.normalized) = rownames(expr)
    ## Rename the normalized expression matrix
    expr = expr.normalized
    
    ## Add the data to the dataObj
    dataObj$expr = expr
    ## Add the genes to the data obj
    dataObj$keys = rownames(expr)
    ## Check that the dataObj is in an acceptable form for metaIntegrator
    stopifnot(checkDataObject(dataObj,"Dataset"))
    ## Add this dataObject to the list of metaObjects
    metaObj[[study]] = dataObj
  }
}

# Wrap the metaObj in a list
metaObj = list(originalData = metaObj)
# Check the the metaObj is in an acceptable form for metaIntegrator
stopifnot(checkDataObject(metaObj, "Meta", "Pre-Analysis"))

# Save the results
saveRDS(metaObj, META.OUT.PATH)