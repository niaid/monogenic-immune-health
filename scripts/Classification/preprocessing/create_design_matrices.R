## Create a list of design matrices and a data frame of corresponding
## sample phenoData in preparation for training the random forest classifiers.
## All matrices should only contain stable features, should be at the subject level,
## and should contain the same subjects as rows in the same order. Features are prepended
## with a name specifying the data type, and rownames are patient ids. 
## The foreground condition group (elsewhere just called the "condition group") 
## and the background condition group are specified via group keys in the file name in the snakemake pipeline.
## Group keys are mapped to vectors of conditions in the script create_classifier_groups.R. 
## Only conditions in the condition group and background condition group are included in the output from this script.

# We load the proper utility functions and packages
source('scripts/util/Processing/averageRepeatSamples.R')
source('scripts/util/Groups/groups.R')
library(Biobase)

# We set the global paths
ESET.IN.PATHS = list(
  somalogic.features = snakemake@input[[1]], #"Data/Somalogic/analysis_output/stability/stable_somalogic_sample_level_features.rds"
  somalogic.modules = snakemake@input[[2]], #"Data/Somalogic/analysis_output/stability/stable_somalogic_sample_level_modules.rds"
  microarray.features = snakemake@input[[3]], #"Data/Microarray/analysis_output/stability/stable_microarray_sample_level_features.rds"
  microarray.modules = snakemake@input[[4]], #"Data/Microarray/analysis_output/stability/stable_microarray_sample_level_modules.rds"
  tbnks = snakemake@input[[5]] #"Data/TBNK/analysis_output/stability/stable_tbnk_sample_level_features.rds"
)

CONDITION.GROUPS.IN.PATH = snakemake@input[[6]] #"Classification/condition_groups.RDS"
BACKGROUND.GROUPS.IN.PATH = snakemake@input[[7]] #"Classification/background_groups.RDS"
SOMALOGIC.MODULES.IN.PATH = snakemake@input[[8]] #"Data/Somalogic/analysis_output/wgcna_results/modules.rds"

DESIGN.MATRIX.OUT.PATH = snakemake@output[[1]] #'Classification/healthy_all_design_matrices_all.RDS'
META.DATA.OUT.PATH = snakemake@output[[2]] #'Classification/healthy_random_forest_sample_meta_data_all.RDS'

# We get the group from the output file name
out.file = basename(DESIGN.MATRIX.OUT.PATH)
out.file = gsub('.RDS$', '', out.file, ignore.case = T)
fields = strsplit(out.file,'_')[[1]]
condition.id = fields[1] # The condition group represents the 'positive' condition for the classifier (e.g. 'healthy', 'cgd', 'xcgd','47cgd','stat1.gof')
background.id = fields[length(fields)] # The background group represents the background pool of all other conditions to consider (e.g. 'PID','AI','all')

# Get the conditons to investigate corresponding to the 'condition' field.
condition.groups = readRDS(CONDITION.GROUPS.IN.PATH)
condition.group = condition.groups[[condition.id]]

# Get the background groups to investiate corresponding to the 'background field'
background.groups = readRDS(BACKGROUND.GROUPS.IN.PATH)
background.group = background.groups[[background.id]]

# We load the relevant data
esets = lapply(ESET.IN.PATHS, readRDS)
somalogic.module.memberships = readRDS(SOMALOGIC.MODULES.IN.PATH)

# We make an eset with just the grey somalogic proteins
somalogic.grey.module = names(somalogic.module.memberships)[somalogic.module.memberships == 'grey']
somalogic.features.eset = esets[['somalogic.features']]
somalogic.grey.eset = somalogic.features.eset[rownames(somalogic.features.eset) %in% somalogic.grey.module, ]
esets[['somalogic.grey']] = somalogic.grey.eset

# We find the visit ids shared by all the data types
patient.ids = Reduce(intersect, lapply(esets, function(x) {x$patient_id}))

# For each eset we
esets = lapply(esets, function(eset) {
  # We subset the data to include only the relevant visits
  eset = eset[, eset$patient_id %in% patient.ids];
  # We average over visit ids from the same patient in the expression sets
  eset = averageRepeatSamples(eset, meta.cols = c('condition','race','gender'))
  # We subset to just the condition of interest and background group
  eset = eset[, eset$condition %in% c(condition.group, background.group)]
  # We rearrange each data set to have subjects in the same order
  eset = eset[, order(eset$patient_id)]
})

# We extract matrices from the esets
Xs = lapply(esets, function(eset) {
  # We create several design matrices with the different features we wish to investigate
  X = t(exprs(eset))
})

# We prefix all of the rownames of each matrix with the revelant data type
Xs = mapply(function(X, name) {colnames(X) = paste(name, colnames(X), sep = '.'); return(X)}, Xs, names(Xs))

# We create a multimodal set with all module scores and tbnks
Xs[['all.modules.with.tbnks']] = cbind(Xs[['microarray.modules']],
                                       Xs[['somalogic.modules']],
                                       Xs[['tbnks']])

# We create a multimodal set with all module scores, grey proteins, and tbnks
Xs[['all.modules.plus.grey.with.tbnks']] = cbind(Xs[['microarray.modules']],
                                                 Xs[['somalogic.modules']],
                                                 Xs[['somalogic.grey']],
                                                 Xs[['tbnks']])

# Here, we add a CBC matrix with just CBC parameters (no lymphocyte phenotyping)
# First, we extract the tbnks matrix
X.tbnks = Xs[['tbnks']]
# We get all the possible absolute and relative features from the lymphocyte populations
tbnk.specific = c('cd3', 'cd4_cd3', 'cd8_cd3', 'cd19', 'nk_cells')
tbnk.specific = c(paste0(tbnk.specific, '_abs'), paste0(tbnk.specific, '_percent'))
tbnk.specific = paste0('tbnks.', tbnk.specific)
# We remove these features from the full set of tbnks features to get cbc-specific features
cbc.specific = setdiff(colnames(X.tbnks), tbnk.specific)
# We subset the tbnks matrix to these features
X.cbcs = X.tbnks[, cbc.specific]
# We add the cbc-specific features matrix to the list of matrices
Xs[['cbcs']] = X.cbcs

# We extract the meta data for the samples
# We choose to use the age from the average of their TBNKs. The ages of a patient
# should be very similar across data types so an arbitrary decision is made.
meta = pData(esets$tbnks)

# We save the data
saveRDS(Xs, file = DESIGN.MATRIX.OUT.PATH)
saveRDS(meta, file = META.DATA.OUT.PATH)
