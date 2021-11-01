## Here we pool together the results of the permutation
## tests to get a final permutation-based pvalue for each feature
## in each classifier.
## The foreground condition group (elsewhere just called the "condition group") 
## and the background condition group are specified via group keys in the file name in the snakemake pipeline.
## Group keys are mapped to vectors of conditions in the script create_classifier_groups.R. 
## Only output from the corresponding RF analysis this condition group and background condition group 
## are included in the output from this script.

## Set globals
## Use the input files from the snakefile to get the permutation results
PERMUTATIONS.IN.PATHS = snakemake@input
PVALS.OUT.PATH = snakemake@output[[1]]#'Classification/results/healthy_rf_pvals_all.RDS'

## Print confirmation that the desired input files are being used
print('Files used for permutation-based pvalue estimation:')
print(PERMUTATIONS.IN.PATHS)

## Load data
results = lapply(PERMUTATIONS.IN.PATHS, readRDS)

## Get the names of the classifiers
classifiers = names(results[[1]])

## For each classifier
p.valss = lapply(classifiers, function(classifier) {
  ## For each permutation test iteration
  p.vals = sapply(results, function(result) {
    ## Get that iteration's pvalues (for each feature)
    ## for that classifier
    result[[classifier]]
  })
  ## Get the average permutation pvalues across
  ## all iterations
  p.vals = rowMeans(p.vals)
  
  ## Here, we set a lower limit for the pvals, reflecting the fact that permutation tests
  ## have precision based on the number of tests run (#tests per iteration x # iterations).
  ## Note that the number of tests per iteration for a classifer equals the number of features
  ## in that classifier, based on how the permutation iterations were run.
  ## For each classifier
  n.tests = length(p.vals) * length(results)
  ## Get the permutation test precision
  precision = 1 / n.tests
  ## For each feature take the larger of the permutation-based pvalue and the precison
  pmax(p.vals, precision)
})

## Name the pvalues by classifier
names(p.valss) = classifiers

## Save results
saveRDS(p.valss, PVALS.OUT.PATH)
