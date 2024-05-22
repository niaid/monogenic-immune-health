## Puts the raw Baltimore Aging Cohort data into an eset

# Load libraries
library(Biobase)

# Set globals
## Baltimore aging cohort hybrid-normalized, calibration-normalized, and median-normalized RFUs
RFUS.IN.PATH = snakemake@input[[1]]#'Reference/ferrucci/raw/CHI-16-018_Hyb.Cal_RFU.txt.adat_RFU.txt'
## Baltimore aging cohort sample metadata
SAMPLES.META.IN.PATH = snakemake@input[[2]]#'Reference/ferrucci/raw/CHI-16-018_Hyb.Cal_RFU.txt.adat_Samples.txt'
## Baltimore aging cohort somamer metadata
SOMAMERS.META.IN.PATH = snakemake@input[[3]]#'Reference/ferrucci/raw/CHI-16-018_Hyb.Cal_RFU.txt.adat_Somamers.txt'

## Baltimore aging cohort eset
ESET.OUT.PATH = snakemake@output[[1]]#'Reference/ferrucci/processed/aging_eset.RDS'

# Load data
RFUs = read.table(RFUS.IN.PATH, sep = '\t', comment.char = '', header = FALSE)
samples.meta = read.table(SAMPLES.META.IN.PATH, sep = '\t', comment.char = '', header = TRUE)
nlines = length(readLines(SOMAMERS.META.IN.PATH)) # There are four lines of comments at the bottom of the file, so we only take the first n - 4 rows
somamers.meta = read.table(SOMAMERS.META.IN.PATH, sep = '\t', comment.char = '',
                           quote = '', header = FALSE, nrow = nlines - 4, row.names = 1)

# Put somamer metadata into a dataframe
somamers.meta = t(somamers.meta)
somamers.meta = as.data.frame(somamers.meta)

# Name the RFUs data frame columns using the somamer ID
# We use somamer IDs rather than the target name to ensure compatibility with the somalogic
# data in the monogenic cohort
colnames(RFUs) = somamers.meta$SomaId

# Name the RFUs data frame forws using the plate ids and positions
ids = paste(samples.meta$PlateId, samples.meta$PlatePosition)
ids = gsub(' ', '_', ids)
rownames(RFUs) = ids

# Convert the RFUs data frame to a matrix
RFUs = as.matrix(RFUs)

# Log transform the data
RFUs = log2(RFUs)

# Create the eset
rownames(somamers.meta) = somamers.meta$SomaId
rownames(samples.meta) = ids
eset = ExpressionSet(t(RFUs))
featureData(eset) = AnnotatedDataFrame(somamers.meta)
phenoData(eset) = AnnotatedDataFrame(samples.meta)

# Save the results
saveRDS(eset, ESET.OUT.PATH)
