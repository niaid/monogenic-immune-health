## Create TBNK training and testing esets with unwanted patients and features removed
## Note that these TBNKs are in the Monogenic metadata database and were pulled from BTRIS.

# Load libraries
library(Biobase)

# Set globals
## Monogenic metadata database
DATABASE.IN.PATH = snakemake@input[[1]]#'Metadata/monogenic.de-identified.metadata.RData'

## the tbnk training eset (prior to cleaning)
TRAINING.OUT.PATH = snakemake@output[[1]]#'Data/TBNK/processed/tbnk_eset_training.rds'
## the tbnk testing eset (prior to cleaning)
TESTING.OUT.PATH = snakemake@output[[2]]#'Data/TBNK/processed/tbnk_eset_testing.rds'

# Load monogenic database and extract tbnks
load(DATABASE.IN.PATH)
tbnks = monogenic.tbnk

# Make Matrix

## Convert patients and visit ids to the 'P_' and 'V_' format
tbnks$patient_id = paste0('P', as.numeric(tbnks$patient_id))
tbnks$visit_id = paste0('V', as.numeric(tbnks$visit_id))

## Make df with only CBC/TBNK parameters (i.e. remove meta data columns)
unit.columns = grep('uom', colnames(tbnks), value = T)
range.columns = grep('range', colnames(tbnks), value = T)
features.to.remove = c(unit.columns, range.columns, c('patient_id', 'visit_id', 'version'))
features = setdiff(colnames(tbnks), features.to.remove)
X = tbnks[, features]

## Name the df rows based on the visit ids
rownames(X) = tbnks$visit_id

## Remove any features that are NA for all samples (i.e. the BANDs)
X = X[, apply(X, 2, function(x) {! all(is.na(x)) })]

## Convert the df to a matrix and transpose it so it is samples x features
X = t(as.matrix(X))

# Make PhenoData

## Choose columns from the monogenic database we wish to use as metadata
meta.features = c('patient_id', 
                  'visit_id',
                  'condition',
                  'patient_age_at_time_of_blood_draw',
                  'gender','race',
                  'ethnicity',
                  'analysis_group',
                  'visit_type')

## Get the corresponding columns from the monogenic database
meta.data = monogenic.all.assays[meta.features]
meta.data = unique(meta.data)
meta.data = meta.data[!is.na(meta.data$visit_id), ]

## Convert patients and visit ids to the 'P_' and 'V_' format
meta.data$visit_id = paste0('V', as.character(meta.data$visit_id))
meta.data$patient_id = paste0('P', as.character(meta.data$patient_id))

## Name the metadata rows based on the visit ids
rownames(meta.data) = meta.data$visit_id

## Select only the visit ids included in the TBNK data
meta.data = meta.data[colnames(X), ]

# Make Feature Data

## Extract units from the units columns in the included features
uoms = sapply(rownames(X), function(feature) {
  name = paste0(feature,'_uom')
  if(name %in% colnames(tbnks)) {
    return(unique(na.omit(tbnks[[name]])))
  } else {
    return(NA)
  }
})

## Make a data frame from the units of measure
f.data = data.frame(uoms = uoms, stringsAsFactors = FALSE)
rownames(f.data) = rownames(X)

## If the units of a feature are missing (i.e. for the TBNKs) and the feature is an absolute quantity,
## make the units that of a concentration
f.data[grepl('_count', rownames(f.data)) & is.na(f.data$uoms), 'uoms'] = '/uL'
## If the units of a feature are missing (i.e. for the TBNKs) and the feature is a relative,
## make the units a percentage
f.data[(!grepl('_count', rownames(f.data))) & is.na(f.data$uoms), 'uoms'] = '%'

# Put together expression set
tbnks.eset = ExpressionSet(X)
phenoData(tbnks.eset) = AnnotatedDataFrame(meta.data)
featureData(tbnks.eset) = AnnotatedDataFrame(f.data)

# Split eset into training and testing
tbnks.train.eset = tbnks.eset[, tbnks.eset$analysis_group == 'Discovery']
tbnks.test.eset = tbnks.eset[, tbnks.eset$analysis_group == 'Validation']

# Save outputs
saveRDS(tbnks.train.eset, file = TRAINING.OUT.PATH)
saveRDS(tbnks.test.eset, file = TESTING.OUT.PATH)
