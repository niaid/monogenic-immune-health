## Create a binary matrix for each visit indicating which medications a patient was on at that visit

# Set globals
## Monogenic metadata database
DATABASE.IN.PATH = snakemake@input[[1]]#'Metadata/monogenic.de-identified.metadata.RData'

## Binary matrix of visits by medications taken, with medications combined into major groups
MEDICATION.TYPES.OUT.PATH = snakemake@output[[1]]#'Medications/medications.types.rds'
## Binary matrix of visits by medications taken, with more resolved medications types
MEDICATION.NAMES.OUT.PATH = snakemake@output[[2]]#'Medications/medications.names.rds'

# Load Database
load(DATABASE.IN.PATH)

# Get data.frame of all medications
medications = monogenic.clinical[c('patient_id', 'visit_id', 'med_type','med_name', 'med_start_date', 'med_end_date')]
medications = unique(medications)

# Change patient_ids and visit_ids to 'P_' and 'V_' format
medications$patient_id = paste0('P', medications$patient_id)
medications$visit_id = paste0('V', medications$visit_id)
monogenic.all.assays$patient_id = paste0('P', monogenic.all.assays$patient_id)
monogenic.all.assays$visit_id = paste0('V', monogenic.all.assays$visit_id)

# Make a map associating visit ids with visit dates
visits = data.frame(visit_date = as.numeric(monogenic.all.assays$visit_date), visit_id = monogenic.all.assays$visit_id)
visits = unique(visits)

visit.dates = visits$visit_date
names(visit.dates) = visits$visit_id

# Associate each visit with its visit date in the medications df 
medications$visit_date = visit.dates[as.character(medications$visit_id)]

# Select only rows for which a medication is being taken during a visit date
select_rows = mapply(function(start, end, date) {
  (is.na(start) & is.na(end)) |
    (is.na(start) & date <= end) | 
    (is.na(end) & start <= date) | 
    (start <= date & date <= end)
}, medications$med_start_date, medications$med_end_date, medications$visit_date)

medications = medications[select_rows, ]

# Associate each visit_id with the corresponding patient id
visit.df = unique(medications[c('patient_id','visit_id')])
visit.map = visit.df$patient_id
names(visit.map) = visit.df$visit_id

# Instantiate a function to turn the long matrix into a binary wide matrix
binarize = function(med_col_name) {
  # Get the visit_ids
  visit_ids = factor(medications$visit_id)
  # Get the medications
  meds = factor(medications[[med_col_name]])
  # Determine if a patient was taking a drug during the visit date
  df.bin = table(visit_id = visit_ids, meds) > 0
  # Convert to a dataframe
  df.bin = as.data.frame(df.bin)
  # Add in the patient_ids and visit_ids
  df.bin$visit_id = rownames(df.bin)
  df.bin$patient_id = visit.map[df.bin$visit_id]
  # Rearrange to make meta data columns first
  df.bin = df.bin[c('patient_id', 'visit_id', levels(meds))]
  # Make all the column names proper names
  colnames(df.bin) = make.names(colnames(df.bin))
  
  return(df.bin)
}

# Do this at both the high and low medications levels
medications.types = binarize('med_type')
medications.names = binarize('med_name')

# Save results
saveRDS(medications.types, MEDICATION.TYPES.OUT.PATH)
saveRDS(medications.names, MEDICATION.NAMES.OUT.PATH)
