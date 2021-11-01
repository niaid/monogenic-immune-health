## This script specifies a map from condition group / background group names to the conditions they encompass. This object is essentially used
## and input to the HI computation pipeline. If you are interested in adding more condition or background groups that the random forest
## can be run on, put them here.

# We load the proper utility functions and packages
source('scripts/util/Groups/groups.R')

# Set globals
CONDITION.GROUPS.OUT.PATH = snakemake@output[[1]]#'Classification/condition_groups.RDS'
BACKGROUND.GROUPS.OUT.PATH = snakemake@output[[2]]#'Classification/background_groups.RDS'

# Create the list of condition groups (the 'positive' class for the random forest)
# The names correspond to ids for the condition groups (abbreviations used in the file names corresponding to the condition groups)
condition.groups = list(
  'cgd' = c('XCGD', '47CGD'),
  'xcgd' = 'XCGD',
  '47cgd' = '47CGD',
  'stat1' = 'STAT1 GOF',
  'job' = 'Job',
  'fmf' = 'FMF',
  'healthy' = 'Healthy'
)

# Create the list of background groups
# It is okay if the background group contains the condition group, they will be eliminated before creating
# the 'negative' class for the classifier
background.groups = list(
  'AI' = util.get_ai(),
  'PID' = util.get_pid(),
  'all' = c(util.get_pid(), util.get_ai(), util.get_tert_terc())
)

# We eliminate NEMO carriers from the background groups
background.groups = lapply(background.groups, function(background.group) {setdiff(background.group, 'NEMO carrier')})

# Save results
saveRDS(condition.groups, CONDITION.GROUPS.OUT.PATH)
saveRDS(background.groups, BACKGROUND.GROUPS.OUT.PATH)
