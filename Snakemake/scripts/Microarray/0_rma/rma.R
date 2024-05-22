library(Biobase)
library(BiocGenerics)
library(pd.hugene.1.0.st.v1)
library(oligo)

#Config ------------------------------------------------
#Input
#Path to input cel files
celDir <- snakemake@input[["cel_dir"]]#"../Monogenic_microarray/Monogenic_Microarray_CEL_File"

#Output
rawOut <- snakemake@output[["raw"]]#"Data/Microarray/raw/raw_featureset.rds"
esetOut <- snakemake@output[["eset"]]#"Data/Microarray/raw/eset_rma_out.rds"

#load data ----------------------------------------------
print("loading Cel files")
celFiles <- list.celfiles(celDir, full.names=TRUE)
rawData <- read.celfiles(celFiles)
print("All files read")

#save raw data ------------------------------------------
saveRDS(rawData, rawOut)

#Perform rma---------------------------------------------
eset <- rma(rawData, target="core")
print("RMA done")

#Save rma expressionset ---------------------------------
saveRDS(eset, esetOut)
print(0)
