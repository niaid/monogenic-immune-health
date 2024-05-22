## Log transform the Somalogic data, put it into an expression set, 
## and export the expression set's training subset, testing subset, and qc subset.
## Note that the Somalogic data has already been normalized by this point via Julian Candia's
## script, but we do not include this in the pipeline for PII reasons.

# load libraries
library(Biobase)

# Set paths
## The hybrid-normalized, calibration-normalized, and median-normalized RFU outputs from the Somalogic assay
SOMALOGIC.RFUS.IN.PATH = snakemake@input[[1]]#'Data/Somalogic/raw/v1/Cal_QC_CHI_Hyb.Cal.MedNorm_RFU.txt'
## Sample metadata associated with the Somalogic assay
SOMALOGIC.SAMPLES.IN.PATH = snakemake@input[[2]]#'Data/Somalogic/raw/v1/Cal_QC_CHI_Samples.txt'
## Somamer metadata associated with the Somalogic assay
SOMALOGIC.SOMAMERS.IN.PATH = snakemake@input[[3]]#'Data/Somalogic/raw/v1/Cal_QC_CHI_Somamers.txt'
## The monogenic metadata database
DATABASE.IN.PATH = snakemake@input[[4]]#'Metadata/monogenic.de-identified.metadata.RData'

## the Somalogic testing eset
TESTING.OUT.PATH = snakemake@output[[1]]#'Data/Somalogic/processed/v1/testing_somalogic.rds'
## the Somalogic training eset
TRAINING.OUT.PATH = snakemake@output[[2]]#'Data/Somalogic/processed/v1/training_somalogic.rds'
## the Somalogic QC eset
QC.OUT.PATH = snakemake@output[[3]]#'Data/Somalogic/processed/v1/qc_somalogic.rds'

# Define the delimiter for reading files
SEP = '\t'

# Load Data
rfus = read.table(SOMALOGIC.RFUS.IN.PATH, sep = SEP, header = FALSE)
sample_metadata = read.table(SOMALOGIC.SAMPLES.IN.PATH, sep = SEP, stringsAsFactors = FALSE, header = TRUE)
somamer_metadata = read.table(SOMALOGIC.SOMAMERS.IN.PATH, sep = SEP, stringsAsFactors = FALSE, header = TRUE)
load(DATABASE.IN.PATH)

# Change data from a dataframe to a matrix
rfus = as.matrix(rfus)

# Transpose data so that rows are features and columns are samples
rfus = t(rfus)

# Add row names to the meta data
rownames(sample_metadata) = paste(sample_metadata$PlateId, sample_metadata$PlatePosition, sep = '-')
rownames(somamer_metadata) = make.names(somamer_metadata$Target)

# Add row names and column names to the Somalogic RFUs
colnames(rfus) = rownames(sample_metadata)
rownames(rfus) = rownames(somamer_metadata)

# Ensure there are no patients in the somalogic data that are not in the database
stopifnot(all(colnames(rfus) %in% rownames(monogenic.somalogic)))

# Get the sample metadata from the Monogenic Database
sample_metadata = monogenic.somalogic[colnames(rfus),]

# Log transform all the somalogic RFUs
log.rfus = log2(rfus)

# Turn into an expression set
somalogic = ExpressionSet(log.rfus)
phenoData(somalogic) = AnnotatedDataFrame(sample_metadata)
featureData(somalogic) = AnnotatedDataFrame(somamer_metadata)

# Add 'V' and 'P' to visit and patient ids respectively
somalogic$patient_id = paste0('P', as.character(somalogic$patient_id))
somalogic$visit_id = paste0('V', as.character(somalogic$visit_id))

# Subset expression set to training cohort
somalogic.train = somalogic[, somalogic$analysis_group == 'Discovery']
somalogic.test = somalogic[, somalogic$analysis_group == 'Validation']
somalogic.qc = somalogic[, somalogic$analysis_group == 'QC']

# Export Training Somalogic
saveRDS(somalogic.train, file = TRAINING.OUT.PATH)
saveRDS(somalogic.test, file = TESTING.OUT.PATH)
saveRDS(somalogic.qc, file = QC.OUT.PATH)