# Load libraries
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(ggpubr)
library(dplyr)
library(tidyr)
library(reshape2)


if(!exists("snakemake")){
  setwd("../../..")
  source("scripts/util/paper/parse_snakemake.R")
  parse_snakemake(rule = "supplemental_figure_4")
}

source("scripts/util/Plotting/tbnk_featurename_replace.R")

# Set paths
## the monogenic metadata database
META.IN.PATH = snakemake@input[[1]]#'Classification/healthy_random_forest_sample_meta_data_all.RDS'
## the LLOO CV results for each patient for each classifier (i.e. the various healthy indexes)
HI.IN.PATH = snakemake@input[[2]]#'Classification/results/healthy_rf_results_all.RDS'
## the results of the JIVE algorithm
JIVE.IN.PATH = snakemake@input[[3]]#"Integration_output/jive/subject/prcomp_list.rds"
## the pvalues for each feature in the random forest classifier
GVI.PVALS.IN.PATH = snakemake@input[[4]]#"Classification/results/healthy_rf_pvals_all.RDS"
## the somalogic subject-level training eset
ESET.IN.PATH = snakemake@input[[5]]#"Data/Somalogic/data_analysis_ready/analysis_ready_subject_level_training_somalogic.rds" 

## The LOO CV results for each patient for each condition-based random forest classifier
HI.CONDITION.IN.PATHS = list(
  CGD = snakemake@input[[6]],#'Classification/results/cgd_rf_results_all.RDS',
  STAT1.GOF = snakemake@input[[7]],#'Classification/results/stat1_rf_results_all.RDS',
  FMF = snakemake@input[[8]],#'Classification/results/fmf_rf_results_all.RDS',
  Job = snakemake@input[[9]]#'Classification/results/job_rf_results_all.RDS'
)

## The meta data used for each condition-based random forest classifier
RF.META.IN.PATHS = list(
  CGD = snakemake@input[[10]],#'Classification/cgd_random_forest_sample_meta_data_all.RDS',
  STAT1.GOF = snakemake@input[[11]],#'Classification/stat1_random_forest_sample_meta_data_all.RDS',
  FMF = snakemake@input[[12]],#'Classification/fmf_random_forest_sample_meta_data_all.RDS',
  Job = snakemake@input[[13]]#'Classification/job_random_forest_sample_meta_data_all.RDS'
)

## The condition-based random forest classifier results for the pvalue of the GVI of each feature via permutation testing
HI.CONDITION.PVALS.IN.PATH = list(
  CGD = snakemake@input[[14]],#'Classification/results/cgd_rf_pvals_all.RDS',
  STAT1.GOF = snakemake@input[[15]],#'Classification/results/stat1_rf_pvals_all.RDS',
  FMF = snakemake@input[[16]],#'Classification/results/fmf_rf_pvals_all.RDS',
  Job = snakemake@input[[17]]#'Classification/results/job_rf_pvals_all.RDS'
)

SUPPLEMENTAL.FIGURE.4a.OUT.PATH = snakemake@output[[1]]#'Paper_1_Figures/Supplemental_Figure_1/S4a.pdf'
SUPPLEMENTAL.FIGURE.4b.OUT.PATH = snakemake@output[[2]]#'Paper_1_Figures/Supplemental_Figure_1/S4b.pdf'
SUPPLEMENTAL.FIGURE.4c.OUT.PATH = snakemake@output[[3]]#'Paper_1_Figures/Supplemental_Figure_1/S4c.pdf'
SUPPLEMENTAL.FIGURE.4d.OUT.PATH = snakemake@output[[4]]#'Paper_1_Figures/Supplemental_Figure_1/S4d.pdf'
SUPPLEMENTAL.FIGURE.4e.OUT.PATH = snakemake@output[[5]]#'Paper_1_Figures/Supplemental_Figure_1/S4e.pdf'
SUPPLEMENTAL.FIGURE.4f.OUT.PATH = snakemake@output[[6]]#'Paper_1_Figures/Supplemental_Figure_1/S4f.pdf'
SUPPLEMENTAL.FIGURE.4g.OUT.PATH = snakemake@output[[7]]#'Paper_1_Figures/Supplemental_Figure_1/S4g.pdf'


# Source utilities
source('scripts/util/Plotting/plot_auc.R')
source('scripts/util/paper/abbrev_cond.R')
# Supplemental figure 4a - Healthy index ROC curves by age group

## Load data
meta = readRDS(META.IN.PATH)
results = readRDS(HI.IN.PATH)

## Extract the HI for each patient
results = results$all.modules.plus.grey.with.tbnks

## Instantiate an empty list to hold ROC data frames
dfs = list()

## Separate the range of ages into three groups
ages = c(0,15,50,100)
groups = c('< 15 yrs','15-50 yrs','> 50 yrs')

## For each age group
for(i in 1:3) {
  
  group = groups[i]
  
  ## Find the subjects from that age group
  select = meta$Age > ages[i] & meta$Age <= ages[i+1]
  
  ## Get the HIs of these subjects
  x = results[select]
  
  ## Get the conditions of these subjects
  y = meta$condition[select]
  
  ## Create an ROC curve for these subjects
  roc = get_roc(x, y, pos = 'Healthy')
  
  ## Get the associated ROC
  auc = get_auc(roc)
  
  ## Store the ROC curve as a data frame, with the age group and auc appended
  df = roc
  df$Age.Group = paste0(group,': ', format(auc, digits = 2))
  
  ## Add the dataframe to the list of data frames
  dfs[[length(dfs) + 1]] = df
}

## Put the dataframes from each age group together
df = Reduce(rbind, dfs)

## Plot the ROC curves
p = ggplot(df, aes(x = fpr, y = tpr, color = Age.Group)) + geom_line() +
  geom_abline(slope = 1, linetype = 'dashed', color = 'grey') + 
  theme_bw() +
  xlab('False Positive Rate') + ylab('True Positive Rate') + labs(color = 'Classifier: AUC') + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = c(.99, .01),
        legend.justification = c(1, 0))

ggsave(SUPPLEMENTAL.FIGURE.4a.OUT.PATH, p, device = 'pdf', height = 5, width = 5)

# Supplemental figure 4b -- ROC curves for the various classifiers

## Load the meta data and random forest LOO CV predictions for each subject
meta = readRDS(META.IN.PATH)
results = readRDS(HI.IN.PATH)

## For each classifier
dfs = mapply(function(x, name) {
  ## Get the ROC curve for that classifier
  roc = get_roc(x = x, y = meta$condition, pos = 'Healthy')
  ## Get the AUC associated with that classifier
  auc = get_auc(roc)
  ## Associate the classifier name and AUC with the ROC dataframe
  roc$Classifier = name
  roc$AUC = auc
  ## Return the ROC dataframe
  return(roc)
}, results, names(results), SIMPLIFY = FALSE)

## Combine the ROC data frames from each classifier
df = Reduce(rbind, dfs)

## Convert the classifier names to be more clear for the plot
conversion = c("microarray.modules" = 'Gene modules', 
               "tbnks" = 'CBC + TBNK',
               "cbcs" = 'CBC',
               "somalogic.modules" = 'Protein modules', 
               "all.modules.with.tbnks" = 'Modules + CBC + TBNK', 
               "all.modules.plus.grey.with.tbnks" = 'Modules + CBC + TBNK + Grey Proteins')
df$Classifier = conversion[df$Classifier]

## Sort the classifiers based on AUC
df = df %>%
  mutate(AUC.formatted = format(df$AUC, digits = 2)) %>%
  mutate(Classifier = paste0(Classifier, ': ', AUC.formatted)) %>%
  arrange(AUC) %>%
  mutate(Classifier = factor(Classifier, levels = unique(Classifier)))

## Plot the various ROC curves together
p = ggplot(df, aes(x = fpr, y = tpr, color = Classifier)) + geom_line() +
  scale_color_manual(values = c('#F8766D', ## we manually choose the colors
                                '#A3A500',
                                '#00BF7D',
                                '#00B0F6',
                                '#E76BF3',
                                '#000000')) +
  geom_abline(slope = 1, linetype = 'dashed', color = 'grey') + ## Include the line y = x to show the performance of a theoretical naive classifier
  theme_bw() +
  xlab('False Positive Rate') + ylab('True Positive Rate') + labs(color = 'Classifier: AUC') + 
  theme(axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.position = c(.99, .01),
        legend.justification = c(1, 0))

ggsave(SUPPLEMENTAL.FIGURE.4b.OUT.PATH, p, device = 'pdf', height = 4, width = 4.5)

# Supplemental Figure 4c -- JIVE versus HI
jive = readRDS(JIVE.IN.PATH) ## Get the jive results
results = readRDS(HI.IN.PATH) ## Get the LOO CV predictions

## Extract the healthy index from the LOO CV predictions
predictions = results[,'all.modules.plus.grey.with.tbnks']
names(predictions) = rownames(results)

## Instantiate a function to create a data frame with healthy index and jive PCs 1:3 from one of the jive matrices
get_df = function(mat.name) {
  pcs = jive[[mat.name]]$x
  ids = intersect(rownames(results), rownames(pcs))
  pcs = pcs[ids, 1:3]
  df = as.data.frame(pcs)
  df$id = rownames(df)
  df = gather(df, key = "variable", value = "value", -id)
  df$predictions = predictions[df$id]
  df$matrix = mat.name
  return(df)
}

mat.names = c('joint','array.ind','soma.ind') ## The various JIVE matrices

## Run the function on each jive matrix
dfs = lapply(mat.names, get_df)
## And stick the results together
df = Reduce(rbind, dfs)

## Order the matrix types by joint, gene, and protein
df <- df %>% 
        mutate(matrix = gsub("joint", "Joint", matrix)) %>%
        mutate(matrix = replace(matrix, matrix == "array.ind", "Transcriptome\nIndividual")) %>%
        mutate(matrix = replace(matrix, matrix == "soma.ind", "Proteome\nIndividual"))

df$matrix = factor(df$matrix, levels = c("Joint", "Transcriptome\nIndividual",
                                         "Proteome\nIndividual"))

## Plot a grid of scatterplots of the correlations between the HI and JIVE PC
## for the joint, protein, and gene JIVE outputs and for PCs 1:3 for each output
p = ggplot(df, aes(x = predictions, y = value)) + ylab('PC Score') +
  xlab('Immune Health Metric') +
  geom_point() + geom_smooth(method = 'lm', formula = y~x, se = FALSE) +
  facet_grid(rows = vars(variable), cols = vars(matrix), scales = "free_y") +
  #xlim(0, .75) +
  stat_cor() +
  theme_bw() + 
  theme(axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        strip.text = element_text(size = 15))

ggsave(SUPPLEMENTAL.FIGURE.4c.OUT.PATH, p, device = 'pdf', height = 5, width = 7)

# Supplemental Figure 4d -- feature-by-feature view of each classifier

## Load the pvalues associated with each feature in the healthy classifiers based on permutation testing of the GVI
results = readRDS(GVI.PVALS.IN.PATH)

## Convert each classifier's name to something more clear (i.e. the data types used by that classifier)
name_conversions = c(cbcs = 'CBC', tbnks = 'CBC + TBNK', microarray.modules = 'Gene Modules', 
                     somalogic.modules = 'Protein Modules', all.modules.with.tbnks = 'Modules + CBC +\nTBNK',
                     all.modules.plus.grey.with.tbnks = 'Modules + CBC +\nTBNK +\nGrey Proteins')

## Instantiate a function to get the negative log10 adjusted pvalue from the permutation test of each feature in the classifier
get.df = function(result, classifier) {
  p.adjusted = p.adjust(result, 'fdr')
  neg.log10.p.adjusted = -1 * log10(p.adjusted)
  df = data.frame(p.adjusted = p.adjusted, neg.log10.pvals = neg.log10.p.adjusted, classifier = classifier)
  df$label = rownames(df)
  return(df)
}

## Run this function on each classifier
dfs = mapply(get.df, results, names(results), SIMPLIFY = FALSE)

## Put all these results together
df = Reduce(rbind, dfs)

## Format the data frames to label significant features
df = df %>%
  mutate(classifier = name_conversions[classifier] %>% factor(levels = name_conversions)) %>%
  mutate(significant = ifelse(p.adjusted < .20, '< .20', '> .20') %>% factor(levels = c('> .20', '< .20'))) %>%
  mutate(label = replace(label, significant == '> .20', ""))

## Manually change some of the labels to be more clear or concise
df = df %>% 
  mutate(label = label %>%
           gsub('tbnks\\.','',.) %>%
           gsub('somalogic\\.grey\\.','',.) %>%
           gsub('somalogic','protein',.) %>%
           gsub('microarray','rna', .) %>%
           gsub("nk_cells_abs", "NK cells(#)", .) %>%
           gsub("nk_cells_percent", "NK cells(%)", .) %>%
           gsub("protein\\.modules\\.purple", "PM2", .) %>%
           gsub('rdw','RDW', .) %>%
           gsub('\\.',' ',.) %>%
           gsub('_',' ',.)
   )

   #%>%
           #gsub('\\.modules','',.) %>%
           #gsub('_',' ',.))

## Instantiate an empty list to hold the plot for each classifier
ps = list()

## Manually choose which color should be used for each classifier (what's science without some art?)
colors = c('CBC' = 'magenta',
           'CBC + TBNK' = 'red',
           'Gene Modules' = 'green',
           'Protein Modules' = 'blue',
           'Modules + CBC +\nTBNK' = 'purple',
           'Modules + CBC +\nTBNK +\nGrey Proteins' = 'black')

## Set the order for the classifier legend
df$classifier = factor(df$classifier, names(colors))

## For each classifier
for(classifier in levels(df$classifier)) {
  ## Get the color for that classifier
  col = colors[[classifier]]
  ## Subset the combined data frame to just the results for that classifier
  df.subset = df %>% filter(classifier == !!classifier)
  ## Create a scatter plot of each classifier (where x is always 0 to keep the points in a vertical line)
  p = ggplot(data = df.subset, aes(x = 0, y = neg.log10.pvals, label = label, color = significant)) + 
    geom_point(show.legend = FALSE) + 
    scale_color_manual(values = c('grey', col)) + # Make the non-significant features grey and the significant ones the chosen color
    xlab(classifier) + 
    ylab('Negative log10 p-value') +
    geom_text_repel(direction='y', nudge_x = .025, hjust = 0, show.legend = FALSE) + # Nudge the feature labels
    xlim(0, .1) + ylim(0, 2.5) + theme_bw()
  
  ## The classifier is the CBCs, which forms the rightmost panel, allow it to display the y-axis title
  if(classifier == 'CBCs') {
    p = p + theme(
      axis.line.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x  = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      axis.line.y = element_line(colour = "black"))
  } else {
    p = p + theme(
      axis.line.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x  = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.title.y = element_blank(),
      axis.line.y = element_line(colour = "black"))
  }
  
  ## Add the plot to the grid
  ps[[length(ps) + 1]] = p
}

## Print the grid
pdf(SUPPLEMENTAL.FIGURE.4d.OUT.PATH, height = 4, width = 11)
grid.arrange(ps[[1]], ps[[2]], ps[[3]], ps[[4]], ps[[5]], ps[[6]], ncol = 6)
dev.off()

# Supplemental Figure 4e -- jPC1-age relationship among top conditions 
jive = readRDS(JIVE.IN.PATH) ## Extract the Jive results
eset = readRDS(ESET.IN.PATH) ## Extract the somalogic eset

## Get the first PC
jPC1 = jive$joint$x[,1]
eset = eset[, names(jPC1)]

## Set the conditions to display
conditions = c('Healthy','47CGD','XCGD','Job','STAT1 GOF','FMF')

## Create a data frame displaying the Age and jPC1 for each condition
df = data.frame(jPC1 = jPC1, condition = eset$condition, age = eset$Age) %>%
  filter(condition %in% conditions) %>%
  mutate(condition = condition %>% as.character %>% abbrev_cond) %>%
  mutate(condition = factor(condition, levels = abbrev_cond(conditions)))

## Create the the Age-PC1 scatterplot for each condition and plot in a grid
p = ggplot(df, aes(x = age, y = jPC1)) + geom_point() + 
  facet_wrap(~condition, ncol = 3, nrow = 2) +
  geom_smooth(method = 'lm', formula = y~x, se = FALSE) +
  ylim(-50, 50) + xlim(0, 80) +
  stat_cor(label.x = 0, label.y = 45) +
  theme_bw() + 
  theme(axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        strip.text = element_text(size = 15))

ggsave(SUPPLEMENTAL.FIGURE.4e.OUT.PATH, p, device = 'pdf', height = 4, width = 6)

# Supplemental Figure 4f -- condition-specific classifiers
results = lapply(HI.CONDITION.IN.PATHS, readRDS) ## Extract the condition-specific classifier results 
metas = lapply(RF.META.IN.PATHS, readRDS) ## Extact the meta data assocaited with each condition-specific classifier

## List the condition groups for each classifier
condition.groups = list(CGD = c('XCGD', '47CGD'),
                        STAT1.GOF = 'STAT1 GOF',
                        FMF = 'FMF',
                        Job = 'Job')

## Create a name conversion map to make the data types underlying each classifier more clear
conversion = c("microarray.modules" = 'Gene modules', 
               "tbnks" = 'CBCs + TBNK',
               "cbcs" = 'CBCs',
               "somalogic.modules" = 'Protein modules', 
               "all.modules.with.tbnks" = 'Modules + CBC + TBNK', 
               "all.modules.plus.grey.with.tbnks" = 'Modules + CBC + TBNK + Grey Proteins')

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
p = ggplot(df, aes(x = variable, y = value, fill = classifier)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  theme_bw() + labs(fill = 'Classifier') + xlab('Condition') + ylab('AUC')

ggsave(SUPPLEMENTAL.FIGURE.4f.OUT.PATH, p, device = 'pdf', height = 6, width = 9)

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

## And plot the associated heatmap using the ggplot tile function
p = ggplot(df, aes(x = x, y = y, fill = NLP)) + geom_tile() + theme_bw() + 
  xlab('Condition') + ylab('Feature') + labs(fill = 'Negative log10 pvalue')
ggsave(SUPPLEMENTAL.FIGURE.4g.OUT.PATH, p, device = 'pdf', height = 6, width = 9)
