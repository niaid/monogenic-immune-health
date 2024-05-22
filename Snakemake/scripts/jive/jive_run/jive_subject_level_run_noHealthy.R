# Perform the JIVE algorithm for all subjects WITH DISEASE (not healthy) 
# with both microarray and somalogic
# and save object with JIVE output matrices for use downstream

library(r.jive)
library(dplyr)
library(Biobase)
library(BiocGenerics)

source("scripts/util/Processing/get_intersecting_data.R")
source("scripts/util/JIVE/JIVE_wrapper.R")

id.col <- "patient_id"

# set paths --------------------------------------------------------- 
#input
ARRAY.ESET.PATH <- snakemake@input[["array_in"]]#"Data/Microarray/analysis_output/stability/stable_microarray_subject_level_features.rds"
SOMA.ESET.PATH <- snakemake@input[["soma_in"]]#"Data/Somalogic/analysis_output/stability/stable_somalogic_subject_level_features.rds"

#output
JIVE.OUT.PATH <- snakemake@output[[1]]#"Integration_output/jive/subject_noHealthy/jive.rds"

if(!dir.exists(dirname(JIVE.OUT.PATH))){
  dir.create(dirname(JIVE.OUT.PATH))
}

# load data ---------------------------------------------------------
array.eset <- readRDS(ARRAY.ESET.PATH)
soma.eset <- readRDS(SOMA.ESET.PATH)

#Subset, remove healthies
array.eset <- array.eset[, array.eset$condition != "Healthy"]
soma.eset <- soma.eset[, soma.eset$condition != "Healthy"]

# get intersecting data ---------------------------------------------
eset.list <- list(array = array.eset, soma = soma.eset)
shared.data <- get_intersecting_data(eset.list) 

#Run Jive ------------------------------------------------------
set.seed(12345) # for reproducibility
jive.out <- jive_wrapper(data.list = shared.data$expr, 
                         z.score = TRUE, frob.scale = TRUE, save.scale.info = TRUE,
                         pdat = shared.data$pdat, 
                         method = "perm",
                         id.col = id.col)



#Save output ----------------------------------------------
saveRDS(jive.out, JIVE.OUT.PATH)
