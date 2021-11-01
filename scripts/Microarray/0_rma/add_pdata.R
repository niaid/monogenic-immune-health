# Takes output from RMA and adds the metadata(phenoData or pData) to an expressionset
# with the normalized expression values.
# subjects no longer in the database are removed.
# patient ID and visit ID are also given P and V prefixes respectively.
# saves expressionset for downstream use.

library(Biobase)
library(BiocGenerics)

#Config ------------------------------------------------
#Input
#Path to input cel files
ESET.IN.PATH <- snakemake@input[["eset"]]#"Data/Microarray/raw/eset_rma_out.rds"
META.DATA.PATH <- snakemake@input[["meta"]]#"Metadata/monogenic.de-identified.metadata.RData"

#Output
esetOut <- snakemake@output[[1]]#"Data/Microarray/raw/eset_rma_with_pData.rds"

#load data ----------------------------------------------
eset <- readRDS(ESET.IN.PATH)
load(META.DATA.PATH)

#Add pData ----------------------------------------------
#There was one patient removed from the protocol that no longer is in monogenic.microarray
eset <- eset[, sampleNames(eset) %in% monogenic.microarray$array_filename]
pdat <- monogenic.microarray[match(sampleNames(eset), monogenic.microarray$array_filename),]
stopifnot(identical(sampleNames(eset), pdat$array_filename))

###Add V and P prefixes to visit ID's
pdat[["patient_id"]] <- paste0("P", pdat[["patient_id"]])
pdat[["visit_id"]] <- paste0("V", pdat[["visit_id"]])

pData(eset) <- pdat

# Save --------------------------------------------------
saveRDS(eset, esetOut)
