# Load libraries
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(ggpubr)
library(dplyr)
library(tidyr)
library(reshape2)

# Source utilities
source('scripts/util/Plotting/plot_auc.R')
source('scripts/util/paper/abbrev_cond.R')

RF.META.IN.PATHS = list(
  CGD = snakemake@input[["cgd_meta"]],#'Classification/cgd_random_forest_sample_meta_data_all.RDS',
  STAT1.GOF = snakemake@input[["stat1_meta"]],#'Classification/stat1_random_forest_sample_meta_data_all.RDS',
  FMF = snakemake@input[["fmf_meta"]],#'Classification/fmf_random_forest_sample_meta_data_all.RDS',
  Job = snakemake@input[["job_meta"]]#'Classification/job_random_forest_sample_meta_data_all.RDS'
)


## The LOO CV results for each patient for each condition-based random forest classifier
HI.CONDITION.IN.PATHS = list(
  CGD = snakemake@input[["cgd_res"]],#'Classification/results/cgd_rf_results_all.RDS',
  STAT1.GOF = snakemake@input[["stat1_res"]],#'Classification/results/stat1_rf_results_all.RDS',
  FMF = snakemake@input[["fmf_res"]],#'Classification/results/fmf_rf_results_all.RDS',
  Job = snakemake@input[["job_res"]]#'Classification/results/job_rf_results_all.RDS'
)


## The condition-based random forest classifier results for the pvalue of the GVI of each feature via permutation testing
HI.CONDITION.PVALS.IN.PATH = list(
  CGD = snakemake@input[["cgd_pvals"]],#'Classification/results/cgd_rf_pvals_all.RDS',
  STAT1.GOF = snakemake@input[["stat1_pvals"]],#'Classification/results/stat1_rf_pvals_all.RDS',
  FMF = snakemake@input[["fmf_pvals"]],#'Classification/results/fmf_rf_pvals_all.RDS',
  Job = snakemake@input[["job_pvals"]]#'Classification/results/job_rf_pvals_all.RDS'
)


## The condition-based random forest classifier results for the pvalue of the GVI of each feature via permutation testing
HI.CONDITION.GVIS.IN.PATH = list(
  CGD = snakemake@input[["cgd_gvis"]],#'Classification/results/cgd_rf_gvis_all.RDS',
  STAT1.GOF = snakemake@input[["stat1_gvis"]],#'Classification/results/stat1_rf_gvis_all.RDS',
  FMF = snakemake@input[["fmf_gvis"]],#'Classification/results/fmf_rf_gvis_all.RDS',
  Job = snakemake@input[["job_gvis"]]#'Classification/results/job_rf_gvis_all.RDS'
)

AUC.FIG.OUT.PATH <- snakemake@output[["auc_fig"]]#"Paper_1_Figures/Supplemental_Figure_2/condition_classifier_auc.pdf"
PVAL.FIG.OUT.PATH <- snakemake@output[["pval_fig"]]

# Supplemental Figure 4f -- condition-specific classifiers
results = lapply(HI.CONDITION.IN.PATHS, readRDS) ## Extract the condition-specific classifier results 

gvis = lapply(HI.CONDITION.GVIS.IN.PATH, readRDS) 

metas = lapply(RF.META.IN.PATHS, readRDS) ## Extact the meta data assocaited with each condition-specific classifier

## List the condition groups for each classifier
condition.groups = list(CGD = c('XCGD', '47CGD'),
                        STAT1.GOF = 'STAT1 GOF',
                        FMF = 'FMF',
                        Job = 'Job')

## Create a name conversion map to make the data types underlying each classifier more clear
conversion = c("microarray.modules" = 'Gene modules', 
               "tbnks" = 'CBCs + Lymphocyte Phenotyping',
               "cbcs" = 'CBCs',
               "somalogic.modules" = 'Protein modules', 
               "all.modules.with.tbnks" = 'Modules + CBCs', 
               "all.modules.plus.grey.with.tbnks" = 'Modules + CBCs + Grey Proteins')

## Insantiate a function to get the AUC associated with each classifier and each condition
get_aucs = function(result, meta, condition.group) {
  ## Get the condition associated with each patient
  conditions = meta[rownames(result), 'condition']
  apply(result, 2, function(x) {
    ## Get the ROC curve associated with each classifier
    roc = get_roc(x = x, y = conditions, pos = condition.group)
    ## Get the AUC of that ROC curve
    get_auc(roc)
  })
}

## Run the function on each of the condition-specific classifier results (and simplify into a matrix)
aucs = mapply(get_aucs, results, metas, condition.groups, SIMPLIFY = T)

## Create a data frame holding the AUCs for each classifier, and melt it
df = as.data.frame(aucs) %>% 
  tibble::rownames_to_column(var = 'classifier') %>%
  mutate(classifier = conversion[classifier]) %>%
  mutate(classifier = factor(classifier, levels = conversion)) %>%
  melt()

## Create grouped barplots for each classifier and each condition
p = df %>% filter(classifier =="Modules + CBCs + Grey Proteins") %>%
        ggplot(aes(x = variable, y = value)) + 
  geom_bar(stat = 'identity') +
  theme_bw() + labs(fill = 'Classifier') + xlab('Condition') + ylab('AUC')

ggsave(AUC.FIG.OUT.PATH, p, device = 'pdf', height = 2, width = 3)


gvis_dat <- lapply(gvis, `[[`, "all.modules.plus.grey.with.tbnks" ) %>%
        lapply(tibble::enframe, name = "feature", value = "gvi") %>%
        bind_rows(.id = "condition") %>%
        mutate(feature = factor(feature), condition = factor(condition))

# Supplemental Figure 4g -- heatmap of gvis
gvi.pvals = lapply(HI.CONDITION.PVALS.IN.PATH, readRDS) ## Get the GVI pvalues associated with each classifier and each condition
pvals = sapply(gvi.pvals, function(x) {x$all.modules.plus.grey.with.tbnks}) ## Extract the pvalues for the features in the classifier with all data types
pvals = as.data.frame(pvals) ## Orangize this matrix into a data frame

## Get the top 5 features from each condition's classifier
top_features = lapply(colnames(pvals), function(group) {
  x = rownames(pvals)[order(pvals[[group]], decreasing = FALSE)]
  x = x[1:5]
})
top_features = unique(unlist(top_features))

## Adjust the pvalues using BH correction within each classifer, and get the negative log 10 adjusted pvalues
pvals = apply(pvals, 2, function(x) {
  x = p.adjust(x, 'fdr')
  x = -log10(x)
})

## Subject to just the top features
pvals = pvals[top_features, ]

## Create an index to associate each feature and each condition with a row and column
n = nrow(pvals)
m = ncol(pvals)
xs = t(matrix(1:m, nrow = m, ncol = n))
ys = matrix(1:n, nrow = m, ncol = n)


## Put the pvalue results, x-indexes, and y-indexes into a data frame
df = data.frame(x = xs[1:(n*m)], y = ys[1:(n*m)], NLP = pvals[1:(n*m)])
df$x = factor(df$x)
levels(df$x) = colnames(pvals)
df$y = factor(df$y)
levels(df$y) = rownames(pvals)

df <- df %>% rename(condition = x, feature = y)

df <- left_join(df, gvis_dat)

df <- df %>%
        mutate(condition = factor(condition, levels = colnames(pvals))) %>%
        mutate(feature = factor(feature, levels = rownames(pvals)))

levels(df$feature) <- levels(df$feature) %>%
        gsub(pattern = "somalogic\\.grey\\.", replacement = "", .) %>%
        gsub(pattern = "microarray\\.modules\\.red", replacement = "TM1 : Inteferon", .)

## And plot the associated heatmap using the ggplot tile function
p = ggplot(df, aes(x = condition, y = feature)) + 
        geom_point(aes(size = NLP, color = gvi)) +
        scale_color_viridis_c() +
        theme_bw() + 
  #xlab('Condition') + ylab('Feature') + 
  labs(size = '-log10(pvalue)') +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1))
ggsave(PVAL.FIG.OUT.PATH, p, device = 'pdf', height = 4, width = 5)
