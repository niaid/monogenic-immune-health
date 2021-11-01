## Formats and reorganizes the CGPs from the OMiCC jamboree,
## combines CGPs from the same study,
## removes unwanted studies,
## and removes unwanted subjects.
## Note that studies actually are considered here to be a 
## GSM-platform pair, so two CGPs from the same GSM but different
## platforms will be left separate. Note that, for error prevention,
## when removing subjects from the CGPs, we do so manually, 
## and check using the series matrices that we removed the correct 
## subjects.

# Set globals
## CGP info
CGPS.IN.PATH = snakemake@input[[1]]#'Reference/jamboree/processed/jamboree_cgps.RDS'
## Series matrix files
SERIES.IN.PATH = snakemake@input[[2]]#'Reference/jamboree/processed/series_matrix_list.rds'

## Cleaned cgps
CGPS.OUT.PATH = snakemake@output[[1]]#'Reference/jamboree/data_analysis_ready/cgps_clean.RDS'

# Read in data
cgps = readRDS(CGPS.IN.PATH)
series = readRDS(SERIES.IN.PATH)

# Remove SLE studies due to the low gene coverage in the OMiCC data
# for this disease
cgps = cgps[names(cgps) != 'SLE']

# Remove unwanted CGPs based on the CGP number for the following reaons:
# DM1_2: Follow up of other study from DM1_1, patients repeated
# RA_3: Largely overlaps with other data already being used in study GSE13501-GPL570
# RA_5: Control group has psoriasis
# RA_10: Largely overlaps with other data already being used in study GSE13501-GPL570
# RA_12: Largely overlaps with other data already being used in study GSE13501-GPL570
# RA_13: Largely overlaps with other data already being used in study GSE13501-GPL570
studies.to.remove = list(DM1 = c(2),
                         MS = c(NULL),
                         RA = c(3, 5, 10, 12, 13),
                         sarcoid = c(NULL))

studies.to.remove = sapply(studies.to.remove, function(study) {paste0('CGP_',study)})

# For each disease
for(disease in names(studies.to.remove)) {
  # For each study in the studies to remove
  studies = studies.to.remove[[disease]]
  cgps.subset = cgps[[disease]]
  for(study in studies) {
    # Remove those studies 
    cgps.subset = cgps.subset[setdiff(names(cgps.subset), studies)]
  }
  # Resave the CGPs for that disease without those studies
  cgps[[disease]] = cgps.subset
}

## Join cgps by study ID, and rename them by study ID

# Copy the CGPs
cgps.new = cgps

# For each disease
for(disease in names(cgps)) {
  # Get the CGPs for that disease
  cgps.subset = cgps[[disease]]
  # Instantiate a list for the study-based groups
  cgps.subset.new = list()
  # For each CGP in the CGPs
  for(cgp in cgps.subset) {
    # Get the study id associated with the cgp
    study.number = cgp$study.info$study
    study.platform = cgp$study.info$platform
    study.id = paste0(study.number,'.',study.platform)
    
    # If the study is not yet in the new cgps list, add it with empty cases and controls
    if(! study.id %in% names(cgps.subset.new)) {
      cgps.subset.new[[length(cgps.subset.new) + 1]] = list(cases = NULL, controls = NULL)
      names(cgps.subset.new)[[length(cgps.subset.new)]] = study.id
    }
    
    # Add in the cases and controls from this cgp to other cases and controls
    # from the same study
    cases = cgps.subset.new[[study.id]]$cases
    cgps.subset.new[[study.id]]$cases = unique(c(cases, cgp$case.names))
    controls = cgps.subset.new[[study.id]]$controls
    cgps.subset.new[[study.id]]$controls = unique(c(controls, cgp$control.names))
  }
  
  # Replace the old OMiCC CGPs with the new ones, organized by study
  cgps.new[[disease]] =  cgps.subset.new
}

# Rename the new cgps
cgps = cgps.new

# Here we edit the studies one at a time to deal with problems that arise

# GSE21942
## Remove technical replicate case samples
gsms.to.remove = c('GSM545843', 'GSM545845')

## We check to make sure we are removing the correct gsms from the series matrix file
# Get the series matrix for this study
check.mat = series$GSE21942.GPL570
# Subset the series matrix to the GSMs have on the study
check.mat = check.mat[check.mat$geo_accession %in% unlist(cgps$MS$GSE21942.GPL570),]
# Get the gsms with 'technical replicate 1' in them in favor of the second tech rep
check.gsms = check.mat$geo_accession[grepl('technical replicate 1', check.mat$title)]
# Check to make sure we removed the correct gsms
stopifnot(all(sort(gsms.to.remove) == sort(check.gsms)))

## We remove these samples
cgps$MS$GSE21942.GPL570$cases = setdiff(cgps$MS$GSE21942.GPL570$cases, gsms.to.remove)
cgps$MS$GSE21942.GPL570$controls = setdiff(cgps$MS$GSE21942.GPL570$controls, gsms.to.remove)

# GSE30210
## We take the last sample from each patient in this longitudinal study
case.gsms.to.keep = c('GSM747681', 'GSM747692', 'GSM747707', 'GSM747725', 'GSM747740', 'GSM747758', 'GSM747766',
                      'GSM747785', 'GSM747800', 'GSM747812', 'GSM747828', 'GSM747841', 'GSM747849', 'GSM747863',
                      'GSM747876', 'GSM747890', 'GSM747903', 'GSM747918')
control.gsms.to.keep = c('GSM747686', 'GSM747695', 'GSM747714', 'GSM747718', 'GSM747732', 'GSM747747', 'GSM747752',
                         'GSM747762', 'GSM747773', 'GSM747793', 'GSM747806', 'GSM747820', 'GSM747835', 'GSM747844',
                         'GSM747854', 'GSM747858','GSM747868', 'GSM747881', 'GSM747899', 'GSM747909', 'GSM747913',
                         'GSM747921')

## We check to make sure we are removing the correct gsms from the series matrix file
# Get the series matrix file
check.mat = series$GSE30210.GPL6947
# Subset to the GSMs we have on the study
check.mat = check.mat[check.mat$geo_accession %in% unlist(cgps$DM1$GSE30210.GPL6947),]
# We remove the time point of the sample name
titles = check.mat$title
titles = sapply(strsplit(titles,'_'), function(x) {x[[1]]})
check.mat$title = titles
# We reverse the order so later samples come first
check.mat = check.mat[rev(1:nrow(check.mat)), ]
# And get the gsms of the non-duplicated names
check.gsms = check.mat$geo_accession[!duplicated(check.mat$title)]
# We check these are equivalent to the gsms we plan to keep
stopifnot(all(sort(c(control.gsms.to.keep, case.gsms.to.keep)) == sort(check.gsms)))

## We keep only these samples
cgps$DM1$GSE30210.GPL6947$cases = case.gsms.to.keep
cgps$DM1$GSE30210.GPL6947$controls = control.gsms.to.keep

# GSE8650
## We remove biological/technical replicates (the last sample is kept for each patient) and patients without symptoms
gsms.to.remove = c('GSM214382', 'GSM214388', 'GSM214390', 'GSM214394', 'GSM214400', 'GSM214406', 'GSM214414',
                   'GSM214416', 'GSM214426', 'GSM214428', 'GSM214442', 'GSM214462', 'GSM214474', 'GSM214484',
                   'GSM214398', 'GSM214420', 'GSM214422', 'GSM214436', 'GSM214438', 'GSM214446', 'GSM214448',
                   'GSM214454', 'GSM214456', 'GSM214458', 'GSM214460', 'GSM214478')

## We check to make sure we are removing the correct gsms from the series matrix file
# We get the series matrix
check.mat = series$GSE8650.GPL96
# We subset it to the GSMs we have from this study
check.mat = check.mat[check.mat$geo_accession %in% unlist(cgps$RA$GSE8650.GPL96),]
# We get the GSMs corresponding to subjects without symptoms
check.gsms.1 = check.mat$geo_accession[check.mat$characteristics_ch1.4 == "Symptoms:  None"]
# We replace the titles with the patient IDs
titles = check.mat$title
titles = sapply(strsplit(titles,'_'), function(x) {x[[3]]})
titles = sapply(strsplit(titles,' '), function(x) {x[[1]]})
titles = tolower(titles)
check.mat$title = titles
# We reverse the order of the matrix to put the last sample for each patient first
check.mat = check.mat[rev(1:nrow(check.mat)), ]
# We get all but the first samples for each patient
check.gsms.2 = check.mat$geo_accession[duplicated(check.mat$title)]
# We put together all the GSMs we've found
check.gsms = unique(c(check.gsms.1, check.gsms.2))
# And we check these are the same as the ones we wish to remove
stopifnot(all(sort(gsms.to.remove) == sort(check.gsms)))

## We remove these samples
cgps$RA$GSE8650.GPL96$cases = setdiff(cgps$RA$GSE8650.GPL96$cases, gsms.to.remove)
cgps$RA$GSE8650.GPL96$controls= setdiff(cgps$RA$GSE8650.GPL96$controls, gsms.to.remove)

## We remove 2 samples that were found to have unreliable diagnoses in the supplemental data from the
## original publication (PMID: 17724127).
## These were JIA patients without diagnosis confirmation upon follow up.
## See supplementary table S1 in original paper for details.
extra.gsms.to.remove = c('GSM214490', 'GSM214492')
cgps$RA$GSE8650.GPL96$cases = setdiff(cgps$RA$GSE8650.GPL96$cases, extra.gsms.to.remove)
cgps$RA$GSE8650.GPL96$controls = setdiff(cgps$RA$GSE8650.GPL96$controls, extra.gsms.to.remove)

# GSE15645
## Remove patients with clinical remission
gsms.to.remove = c('GSM391602', 'GSM391603', 'GSM391604', 'GSM391605', 'GSM391606',
                   'GSM391607', 'GSM391608', 'GSM391609', 'GSM391610', 'GSM391611',
                   'GSM391612', 'GSM391613', 'GSM391614', 'GSM391615', 'GSM391616')

## Check to make sure we are removing the correct patients
# Get series matrix
check.mat = series$GSE15645.GPL570
# Subset to the gsms we have
check.mat = check.mat[check.mat$geo_accession %in% unlist(cgps$RA$GSE15645.GPL570),]
# Find all samples with CR in the title
titles = check.mat$title
# 'CR' stands for clinical remission and 'CRM' stands for clinical remission with medication
check.gsms = check.mat$geo_accession[grepl('_CR_', titles) | grepl('_CRM_', titles)]

# Ensure we are removing the right samples
stopifnot(all(sort(gsms.to.remove) == sort(check.gsms)))

## We remove these patients
cgps$RA$GSE15645.GPL570$cases = setdiff(cgps$RA$GSE15645.GPL570$cases, gsms.to.remove)
cgps$RA$GSE15645.GPL570$controls = setdiff(cgps$RA$GSE15645.GPL570$controls, gsms.to.remove)

# GSE18781
## Case and control GSMs were flipped here (in OMiCC)
controls = cgps$sarcoid$GSE18781.GPL570$cases
cases = cgps$sarcoid$GSE18781.GPL570$controls

cgps$sarcoid$GSE18781.GPL570$cases = cases
cgps$sarcoid$GSE18781.GPL570$controls = controls

# GSE42834
## Remove patients with non-active sarcoid
gsms.to.remove = c('GSM1050754', 'GSM1050759', 'GSM1050762', 'GSM1050763', 'GSM1050766', 'GSM1050774',
                   'GSM1050780', 'GSM1050783', 'GSM1050789', 'GSM1050793', 'GSM1050797', 'GSM1050816',
                   'GSM1050843', 'GSM1050864', 'GSM1050931', 'GSM1050933', 'GSM1050949', 'GSM1050969',
                   'GSM1050973', 'GSM1050975', 'GSM1050976', 'GSM1050977')

## Check to make sure we are removing the correct patients
check.mat = series$GSE42834.GPL10558
check.mat = check.mat[check.mat$geo_accession %in% unlist(cgps$sarcoid$GSE42834.GPL10558), ]
check.gsms = check.mat[check.mat$characteristics_ch1.2 == 'disease state: Non-active sarcoidosis', 'geo_accession']
stopifnot(sort(check.gsms) == sort(gsms.to.remove))

## Remove these patients
cgps$sarcoid$GSE42834.GPL10558$cases = setdiff(cgps$sarcoid$GSE42834.GPL10558$cases, gsms.to.remove)
cgps$sarcoid$GSE42834.GPL10558$controls = setdiff(cgps$sarcoid$GSE42834.GPL10558$controls, gsms.to.remove)

# Check that no GSMs are repeated
gsms = unname(unlist(cgps, recursive = T))
stopifnot(max(table(gsms)) == 1)

# Save the results
saveRDS(cgps, CGPS.OUT.PATH)
