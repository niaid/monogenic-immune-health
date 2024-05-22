## Clean and format training TBNKs, selecting the desired samples and features, 
## adding in the neutrophil lymphocyte ratio (NLR) manually. Create subject and sample level esets.
## Convert all features expressed in percent lymphocytes to percent of WBCs (and also rename them with 'percent')
## Manually remove some additional features for the following reasons: 
## 1) immature_granulocytes_abs: many samples were NA
## 2) immature_granulocytes_percent: many samples were NA
## 3) mpv: many samples were NA
## 4) nucleated_rbc_abs: too few counts to be reliable
## 5) nucleated_rbc_percent: too few counts to be reliable 

# Load relevant libraries
library(Biobase)

# Source utilities
source('scripts/util/Processing/averageRepeatSamples.R')

# Set globals
## the tbnk training eset (prior to cleaning)
TBNKS.IN.PATH = snakemake@input[[1]]#'Data/TBNK/processed/tbnk_eset_training.rds'

## the cleaned up sample-level tbnk training eset
TBNKS.SAMPLE.LEVEL.OUT.PATH = snakemake@output[[1]]#'Data/TBNK/data_analysis_ready/tbnk_training_sample_level_eset.rds'
## the cleaned up subject-level tbnk training eset
TBNKS.SUBJECT.LEVEL.OUT.PATH = snakemake@output[[2]]#'Data/TBNK/data_analysis_ready/tbnk_training_subject_level_eset.rds'

# Load data
tbnks = readRDS(TBNKS.IN.PATH)

# Remove unwanted samples for the following reasons
# V316 & V318: Lymphocyte subpopulations reported in TBNKs do not sum to total lymphocytes
# V40: Lymphocyte percentages reported in TBNKs are not equal to absolute population counts divided by total lymphocytes
outliers = c('V40', 'V316','V318')
tbnks = tbnks[, ! tbnks$visit_id %in% outliers]

# Convert all features expressed in percent lymphocytes
# to percent of WBCs (and also rename them with 'percent')
lymphocyte.subsets = c('cd3', 'cd4_cd3', 'cd8_cd3', 'cd19', 'nk_cells')
lymphocytes.percent = exprs(tbnks)['lymphocytes_percent', ]
for(feature in lymphocyte.subsets) {
  exprs(tbnks)[feature, ] = exprs(tbnks)[feature, ] * lymphocytes.percent / 100
  rownames(tbnks)[rownames(tbnks) == feature] = paste0(feature, '_percent')
}

# Rename _count to _abs for lymphocytes
rownames(tbnks) = gsub('_count', '_abs', rownames(tbnks))

# Choose the features we wish to include (note platelets are included as a 'population')
features = rownames(tbnks)
marrow.features = c('wbc', 'rbc', 'hemoglobin', 'mcv', 'mch', 'mchc', 'rdw')
absolute.population.features = grep('_abs', features, value = T)
relative.population.features = grep('_percent', features, value = T)
features.to.include = c(marrow.features, absolute.population.features, relative.population.features)

# Manually remove some additional features for the following reasons
# immature_granulocytes_abs: many samples were NA
# immature_granulocytes_percent: many samples were NA
# mpv: many samples were NA
# nucleated_rbc_abs: too few counts to be reliable
# nucleated_rbc_percent: too few counts to be reliable 
features.to.exclude = c('immature_granulocytes_abs',
                        'immature_granulocytes_percent',
                        'nucleated_rbc_abs',
                        'nucleated_rbc_percent',
                        'mpv')
features.to.include = setdiff(features.to.include, features.to.exclude)

# Subset to the desired features
tbnks = tbnks[features.to.include,]

# Only use samples in which all desired features have complete data
tbnks = tbnks[, complete.cases(t(exprs(tbnks)))]

# Add in the neutrophil to lymphocyte ratio to the features and feature data
df = as.data.frame(t(exprs(tbnks)))
df$'NLR' = df[['neutrophil_abs']] / df[['lymphocytes_abs']] * 100
f.data = fData(tbnks)
f.data['NLR',] = '%'

# Convert the data frame back to a matrix and transpose
X = t(as.matrix(df))

# Wrap everything up into an expression set
tbnks.samples = ExpressionSet(X)
phenoData(tbnks.samples) = phenoData(tbnks)
featureData(tbnks.samples) = AnnotatedDataFrame(f.data)

# Average over biological repeats
tbnks.subjects = averageRepeatSamples(tbnks.samples)

# Save esets
saveRDS(tbnks.samples, TBNKS.SAMPLE.LEVEL.OUT.PATH)
saveRDS(tbnks.subjects, TBNKS.SUBJECT.LEVEL.OUT.PATH)
