## Gets DE statistics using the limma model fits. Statistics included are t-statistic, linear model effect size, cohen's D statistic
## p value, adjusted pvalue (FDR), and stable feature adjust pvalue (FDR corrected only among the stable features)

# Set seed
set.seed(1704)

# Load libraries and source utilities
library(limma)
library(Biobase)
library(dplyr)
source('scripts/util/DifferentialExpression/limma.R')

# Set globals
## Variance partitions for accessing stable features
VP.IN.PATHS = list(
  somalogic.modules = snakemake@input[[1]],#'Data/Somalogic/analysis_output/stability/somalogic_modules_standard_vp.rds',
  somalogic.features = snakemake@input[[2]],#'Data/Somalogic/analysis_output/stability/somalogic_features_standard_vp.rds',
  microarray.modules = snakemake@input[[3]],#'Data/Microarray/analysis_output/stability/microarray_modules_standard_vp.rds',
  microarray.features = snakemake@input[[4]],#'Data/Microarray/analysis_output/stability/microarray_features_standard_vp.rds',
  tbnks.features = snakemake@input[[5]]#'Data/TBNK/analysis_output/stability/tbnks_features_standard_vp.rds'
)

## DE model fits
FIT.IN.PATHS = list(
  somalogic.modules = snakemake@input[[6]],#'Data/Somalogic/analysis_output/differential_expression/somalogic_modules_DE_fit.rds',
  somalogic.features = snakemake@input[[7]],#'Data/Somalogic/analysis_output/differential_expression/somalogic_features_DE_fit.rds',
  microarray.modules = snakemake@input[[8]],#'Data/Microarray/analysis_output/differential_expression/microarray_modules_DE_fit.rds',
  microarray.features = snakemake@input[[9]],#'Data/Microarray/analysis_output/differential_expression/microarray_features_DE_fit.rds',
  tbnks.features = snakemake@input[[10]]#'Data/TBNK/analysis_output/differential_expression/tbnks_features_DE_fit.rds'
)

## Where to place the statistics derived in this script (just stable features)
RESULT.OUT.PATHS = list(
  somalogic.modules = snakemake@output[[1]],#'Data/Somalogic/analysis_output/differential_expression/somalogic_modules_DE_results.rds',
  somalogic.features = snakemake@output[[2]],#'Data/Somalogic/analysis_output/differential_expression/somalogic_features_DE_results.rds',
  microarray.modules = snakemake@output[[3]],#'Data/Microarray/analysis_output/differential_expression/microarray_modules_DE_results.rds',
  microarray.features = snakemake@output[[4]],#'Data/Microarray/analysis_output/differential_expression/microarray_features_DE_results.rds',
  tbnks.features = snakemake@output[[5]]#'Data/TBNK/analysis_output/differential_expression/tbnks_features_DE_results.rds'
)

# Where to place the intermediate results (with all features)
INTERMEDIATE.OUT.PATHS = list(
  somalogic.modules = snakemake@output[[6]],#'Data/Somalogic/analysis_output/differential_expression/somalogic_modules_DE_intermediates.rds',
  somalogic.features = snakemake@output[[7]],#'Data/Somalogic/analysis_output/differential_expression/somalogic_features_DE_intermediates.rds',
  microarray.modules = snakemake@output[[8]],#'Data/Microarray/analysis_output/differential_expression/microarray_modules_DE_intermediates.rds',
  microarray.features = snakemake@output[[9]],#'Data/Microarray/analysis_output/differential_expression/microarray_features_DE_intermediates.rds',
  tbnks.features = snakemake@output[[10]]#'Data/TBNK/analysis_output/differential_expression/tbnks_features_DE_intermediates.rds'
)

# Load fits
fits = lapply(FIT.IN.PATHS, readRDS)

# Load vps
vps = lapply(VP.IN.PATHS, readRDS)

# Get list of condition groups we wish to analyze (note that the *CONDITIONS variables are globals
# from the utility script)
groups = c(make_condition_groups(conditions = CONDITIONS),
           make_AI_PID_groups(ai.conditions = AI_CONDITIONS,
                              pid.conditions = PID_CONDITIONS,
                              tert.terc.conditions = TERT_TERC_CONDITIONS))

# Instantiate a function to get stats associated with stable features of an eset
run_stats = function(fit, vp) {
  
  # For the versus-healthy and versus-all analysis
  stats = lapply(c(FALSE, TRUE), function(cross) {
    
    # Design contrast matrix using the specified groups
    contrast = make_contrasts_mat(fit, groups, cross)
    
    # Fit the contrast matrix
    contrast.fit = contrasts.fit(fit, contrast)
    
    # Use ebayes or traditional t-test based on number of features
    if(nrow(vp) < 100) {
      stats = get_traditional_stats(contrast.fit)
    } else {
      stats = get_ebayes_stats(contrast.fit)
    }
    
    # Add an adjusted pvalue using FDR
    stats$adj.P.Val = apply(stats$P.Value, 2, function(x) {
      p.adjust(x, 'fdr')
    })
    
    return(stats)
    
  })
  
  names(stats) = c('versus.healthy', 'versus.all')
  return(stats)
}

# Instantiate a function to subset statistics from a data type to just the stable features
subset_stats = function(statss, vp) {
  
  # Get stable feautres
  stable.features = vp@row.names[vp$Patient >= .5]
  
  # For each of the options (versus healthy and versus all)
  statss.stable = lapply(statss, function(stats) {
    
    # For each statistic
    stats.stable = lapply(stats, function(stat) {
      
      # Subset the statistic to the stable features
      stat[stable.features, ]
      
    })
    
    # Add an adjusted pvalue, only among the stable features
    stats.stable$stable.feature.adj.P.Val = apply(stats.stable$P.Value, 2, function(x) {
      
      p.adjust(x, 'fdr')
      
    })
    
    
    return(stats.stable)
    
  })
}

# Apply the stat-computing function over all esets and vps
results = mapply(run_stats, fits, vps, SIMPLIFY = FALSE)

# Apply the stability filtering function over all esets and vps
results.stable = mapply(subset_stats, results, vps, SIMPLIFY = FALSE)

# Save intermediates
mapply(function(result, out.path) {
  saveRDS(result, out.path)
}, results, INTERMEDIATE.OUT.PATHS)

# Save stats
mapply(function(result, out.path) {
  saveRDS(result, out.path)
}, results.stable, RESULT.OUT.PATHS)