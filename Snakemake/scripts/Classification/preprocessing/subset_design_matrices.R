## Here we create a list of design matrices we wish to use in the random forests,
## eliminating options like the microarray features as they are too large.
## The foreground condition group (elsewhere just called the "condition group") 
## and the background condition group are specified via group keys in the file name in the snakemake pipeline.
## Group keys are mapped to vectors of conditions in the script create_classifier_groups.R. 
## Only conditions in the condition group and background condition group are included in the output from this script.

# Set paths
DESIGN.MATRICES.IN.PATH = snakemake@input[[1]]#'Classification/healthy_design_matrices_all.RDS'
META.DATA.IN.PATH = snakemake@input[[2]]#'Classification/healthy_random_forest_sample_meta_data_all.RDS'

DESIGN.MATRICES.OUT.PATH = snakemake@output[[1]]#'Classification/healthy_random_forest_design_matrices_all.RDS'

# Read in design matrices and meta data
Xs = readRDS(DESIGN.MATRICES.IN.PATH)
meta = readRDS(META.DATA.IN.PATH)

# Subset the design matrices to those we wish to use for the random forest
features = c('cbcs',
             'tbnks',
             'microarray.modules',
             'somalogic.modules',
             'all.modules.with.tbnks',
             'all.modules.plus.grey.with.tbnks')

Xs = Xs[features]

# Save design matrices and meta data
saveRDS(Xs, DESIGN.MATRICES.OUT.PATH)