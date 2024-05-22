# This script computes the PCA for a given jive run and saves them as a list. 
# The PCA is computed for joint as well as the individual matrices.
# The output list contains prcomp objects for the joint and the two individual matrices
# It uses file path objects created by snakemake in order to access the proper in/output paths

#Source Scripts ----------------------------------------------
source("scripts/util/JIVE/JIVE_pca.R")

#Inputs -----------------------------------------------------
JIVE.PATH <- snakemake@input[["jive_obj"]]

#Outputs ----------------------------------------------------
OUT.PRCOMP.LIST.PATH <- snakemake@output[[1]]

#load data ---------------------------------------------------------
jive <- readRDS(JIVE.PATH)

#Run PCA ------------------------------------------------------------------
prcomp.list <- get_jive_pca(jive)

#Flip PC1 of the joint so that correlated with healthy index
prcomp.list$joint$x[, "PC1"] <- -prcomp.list$joint$x[, "PC1"]

#Save --------------------------------------------------------------
saveRDS(prcomp.list, OUT.PRCOMP.LIST.PATH)
