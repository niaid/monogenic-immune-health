# Load packages
library(ggplot2)
library(dplyr)
library(tidyr)
source('scripts/util/paper/abbrev_cond.R')
source('scripts/util/Groups/groups.R')
source('scripts/util/Plotting/plot_auc.R')

# Set paths
## The healthy index using AI as a background
AI.BASED.HI.IN.PATH = snakemake@input[[1]] #'Classification/results/healthy_rf_results_AI.RDS'
## The healthy index feature GVI pvalues from permutation testing using AI as a background
AI.BASED.HI.PVALS.IN.PATH = snakemake@input[[2]] #'Classification/results/healthy_rf_pvals_AI.RDS'
## The AI and healthy sample meta data
AI.BASED.HI.META.IN.PATH = snakemake@input[[3]] #'Classification/healthy_random_forest_sample_meta_data_AI.RDS'
## The healthy index using PID as a background
PID.BASED.HI.IN.PATH = snakemake@input[[4]] #'Classification/results/healthy_rf_results_PID.RDS'
## The healthy index feature GVI pvalues from permutation testing using PID as a background
PID.BASED.HI.META.IN.PATH = snakemake@input[[5]] #'Classification/healthy_random_forest_sample_meta_data_PID.RDS'
## The PID and healthy sample meta data
PID.BASED.HI.PVALS.IN.PATH = snakemake@input[[6]] #'Classification/results/healthy_rf_pvals_PID.RDS'
## The healthy index using all subjects as a background
ALL.BASED.HI.IN.PATH = snakemake@input[[7]] #'Classification/results/healthy_rf_results_all.RDS'
## The healthy index feature GVI pvalues from permutation testing using all conditions as a background
ALL.BASED.HI.META.IN.PATH = snakemake@input[[8]] #'Classification/healthy_random_forest_sample_meta_data_all.RDS'
## The all conditiion and healthy sample meta data
ALL.BASED.HI.PVALS.IN.PATH = snakemake@input[[9]] #'Classification/results/healthy_rf_pvals_all.RDS'

## The PID predictions after training the classifier on just AI and Healthy
PID.PREDICTIONS.FROM.AI.IN.PATH = snakemake@input[[10]] #'Classification/predictions/healthy_rf_PID_predictions_using_AI_index.RDS'
## The AI predictions after training the classifier on just PID and Healthy
AI.PREDICTIONS.FROM.PID.IN.PATH = snakemake@input[[11]] #'Classification/predictions/healthy_rf_AI_predictions_using_PID_index.RDS'

## The figure out path
PDF.OUT.PATH = snakemake@output[[1]] #'Paper_1_Figures/Supplemental_Figure_4/figure_4_AI_and_PID_HI_addendum.pdf'

# Get group constituents
AI = util.get_ai()
PID = util.get_pid()
Telo = util.get_tert_terc()

# Instantiate a function for creating AUC plots
make_auc_plot = function(results, meta, title) {
  result = results$all.modules.plus.grey.with.tbnks
  roc = get_roc(result, meta$condition, 'Healthy')
  auc = get_auc(roc)
  
  auc = format(auc, digits = 2)
  
  p = ggplot(roc, aes(x = fpr, y = tpr)) + geom_line(color = 'black', show.legend = FALSE) + 
    geom_abline(slope = 1, linetype = 'dashed', color = 'grey') + 
    theme_bw() + geom_text(aes(x = .75, y = .25), size = 4, label = paste0('AUC: ', auc), show.legend = FALSE) + 
    xlab('False Positive Rate') + ylab('True Positive Rate') + 
    ggtitle(title) +
    theme(axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8))
}

# Instantiate a function for making the bar plots
make_bar_plots = function(results, meta, title) {
  df = data.frame(healthy.index = results$all.modules.plus.grey.with.tbnks, 
                  condition = meta$condition) %>%
    mutate(condition = as.character(condition)) %>%
    mutate(group = condition %>% # Add in condition super-type
             replace(condition %in% AI, 'AI') %>%
             replace(condition %in% PID, 'PID') %>%
             replace(condition %in% Telo, 'Telo')) %>%
    mutate(group = factor(group, levels = intersect(c('Healthy','AI','Telo','PID'), unique(group)))) %>%
    mutate(condition = abbrev_cond(condition)) # We use the abbreviated condition names
  
  ## Compute the median healthy index for each condition
  condition.median.healthy.indexes = df %>%
    group_by(condition) %>%
    summarise(condition.median.healthy.index = median(healthy.index))
  
  ## Add in the median healthy index for each condition to the original data frame
  df = df %>%
    right_join(condition.median.healthy.indexes, by = 'condition') %>%
    arrange(as.numeric(group), desc(condition.median.healthy.index)) %>% # Sort by condition super-type and then median healthy index
    mutate(condition = factor(condition, condition %>% unique)) %>%
    mutate(condition = relevel(condition, abbrev_cond('Healthy'))) %>% # Make sure Healthy is the first level
    mutate(condition = factor(condition, levels = rev(levels(condition))))
  
  ## Plot the bar plots
  HI_max = max(df$healthy.index) + .01
  p1 = ggplot(df, aes(x = condition, y = healthy.index, fill = group)) + 
    geom_boxplot(outlier.colour = NA) + 
    ylim(0, HI_max) +
    theme_bw() + geom_jitter() + coord_flip() + ggtitle(title) +
    theme(axis.text.x = element_text(size = 15),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.key.size = unit(2,"line"))
}

# Instantiate a function for plotting the top features from each classifier
make_pval_plots = function(p.vals, title) {
  p.vals = p.vals$all.modules.plus.grey.with.tbnks
  p.adjusted = p.adjust(p.vals, 'fdr')
  neg.log10.p.adjusted = -1 * log10(p.adjusted)
  
  ## Create a data frame with the feature names, data type, and negative log 10 pvalues, just for features passing the FDR cutoff
  df = data.frame(label = names(p.vals), p.adjusted = p.adjusted, neg.log10.pvals = neg.log10.p.adjusted) %>%
    filter(p.adjusted < .20) %>%
    select(-p.adjusted) %>%
    mutate(data.type = label %>%
             as.character() %>%
             replace(., grepl('somalogic\\.grey\\.', .), 'Grey\nModule\nProteins') %>%
             replace(., grepl('somalogic\\.modules\\.', .), 'Protein\nModule\nScores') %>%
             replace(., grepl('tbnks\\.', .), 'CBC +\nLymphocyte\nPhenotyping')) %>%
    mutate(label = label %>%
             as.character() %>%
             gsub('somalogic\\.grey\\.', '', .) %>%
             gsub('somalogic\\.modules\\.', '', .) %>%
             gsub('microarray\\.modules\\.', '', .) %>%
             gsub('tbnks\\.', '', .))
  
  ## We now manually clean up the feature names one-by-one to make them look better when plotting
  df = df %>%
    mutate(label = label %>%
             gsub('nk_cells_percent','NK Cells %', .) %>%
             gsub('nk_cells_abs','# NK Cells', .) %>%
             gsub('MIP.1a','MIP 1a', .) %>%
             gsub('purple','Proteomic Purple Module', .) %>%
             gsub('Cathepsin.H','Cathepsin H', .) %>%
             gsub('IL.18.Ra','IL-18 Receptor 1', .) %>%
             gsub('rdw','RDW', .) %>%
             gsub('LD78.beta','LD78 beta', .))
  
  ## We order the features by negative log 10 pvalue
  df = df %>% 
    arrange(neg.log10.pvals) %>%
    mutate(label = factor(label, levels = label)) %>%
    mutate(data.type = factor(data.type))
  
  ## We plot the feature p-values in a bar plot
  p = ggplot(df, aes(y = neg.log10.pvals, x = label, fill = data.type)) + 
    geom_bar(stat="identity") + theme_bw() + xlab('Parameter') + ylab('Negative log10 q-values') + coord_flip() +
    scale_fill_manual(values = c('darkblue','steelblue', 'lightblue')) + labs(fill = 'Data Type') +
    geom_hline(aes(yintercept = -log10(.20), linetype = 'FDR = .20'), color = 'gray', size = 1) +
    scale_linetype_manual(values = 'dashed') + ggtitle(title) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10))
}

# Instantiate a function to plot correlations between healthy indexes in all subjects
make_scatter_plots_one_group = function(results.1, results.2, 
                                        meta.1, meta.2, 
                                        title,
                                        x.label, y.label) {
  
  overlapping.subjects = intersect(rownames(results.1), rownames(results.2))
  
  results.1 = results.1[overlapping.subjects, ]
  meta.1 = meta.1[overlapping.subjects, ]
  results.2 = results.2[overlapping.subjects, ]
  meta.2 = meta.2[overlapping.subjects, ]
  
  stopifnot(as.character(meta.1$condition) == as.character(meta.2$condition))
  meta = meta.1
  
  hi.df = data.frame(HI.1 = results.1$all.modules.plus.grey.with.tbnks,
                     HI.2 = results.2$all.modules.plus.grey.with.tbnks,
                     condition = meta$condition)
  
  p = ggplot(hi.df, aes(x = HI.1, y = HI.2)) + geom_smooth(method = 'lm', formula = y ~ x, se = F) + 
    geom_point(aes(color = condition)) + ggpubr::stat_cor() +
    theme_bw() +
    ggtitle(title) + xlab(x.label) + ylab(y.label)
}

# Instantiate a function to plot correlations between healthy indexes in cases and controls
make_scatter_plots_two_group = function(results.1, results.2, 
                                        meta.1, meta.2, 
                                        title,
                                        group1 = 'Healthy', 
                                        group1.name = 'Healthy Control',
                                        group2.name,
                                        x.label, y.label) {
  
  overlapping.subjects = intersect(rownames(results.1), rownames(results.2))
  
  results.1 = results.1[overlapping.subjects, ]
  meta.1 = meta.1[overlapping.subjects, ]
  results.2 = results.2[overlapping.subjects, ]
  meta.2 = meta.2[overlapping.subjects, ]
  
  stopifnot(as.character(meta.1$condition) == as.character(meta.2$condition))
  meta = meta.1
  
  hi.df = data.frame(HI.1 = results.1$all.modules.plus.grey.with.tbnks,
                     HI.2 = results.2$all.modules.plus.grey.with.tbnks,
                     condition = meta$condition)
  hi.df = hi.df %>% 
    mutate(
      group = ifelse(hi.df$condition %in% group1, group1.name, group2.name) %>%
        factor(levels = c(group1.name, group2.name))
    )
  
  p = ggplot(hi.df, aes(x = HI.1, y = HI.2, group = group)) + geom_smooth(method = 'lm', formula = y ~ x, se=F) + 
    geom_point(aes(color = group)) + ggpubr::stat_cor(aes(group = group, color = group)) +
    theme_bw() +
    ggtitle(title) + xlab(x.label) + ylab(y.label)
}

# Instantiate a function to plot correlations between healthy indexes in a new group
make_scatter_plots_multi_group = function(results.1, results.2, 
                                  meta.1, meta.2, 
                                  title, x.label, y.label, remove.healthy = T) {
  
  overlapping.subjects = intersect(rownames(results.1), rownames(results.2))
  
  results.1 = results.1[overlapping.subjects, ]
  meta.1 = meta.1[overlapping.subjects, ]
  results.2 = results.2[overlapping.subjects, ]
  meta.2 = meta.2[overlapping.subjects, ]
  
  stopifnot(as.character(meta.1$condition) == as.character(meta.2$condition))
  meta = meta.1
  
  hi.df = data.frame(HI.1 = results.1$all.modules.plus.grey.with.tbnks,
                     HI.2 = results.2$all.modules.plus.grey.with.tbnks,
                     condition = meta$condition)
  
  if(remove.healthy) {
    hi.df = hi.df %>% filter(condition != "Healthy")
  }
  
  hi.df = hi.df %>% 
    filter(condition %in% names(table(condition))[table(condition) > 5])
  p = ggplot(hi.df, aes(x = HI.1, y = HI.2, group = condition)) + 
    geom_smooth(aes(color = condition), method = 'lm', formula = y ~ x, se = F) + 
    geom_point(aes(color = condition)) + ggpubr::stat_cor(aes(group = condition, color = condition)) +
    theme_bw() +
    ggtitle(title) + xlab(x.label) + ylab(y.label)
}

# Start pdf
pdf(PDF.OUT.PATH)

# Addendum Figure 1 -- AUC performance of the original healthy index among AI patients and controls
results = readRDS(ALL.BASED.HI.IN.PATH)
meta = readRDS(ALL.BASED.HI.META.IN.PATH)
results = results[meta$condition %in% c(AI, 'Healthy'), ]
meta = meta[meta$condition %in% c(AI, 'Healthy'), ]
title = 'ROC Curve of Original Classifier LOO CV predictions\namong AI and Healthy Subjects'

p = make_auc_plot(results, meta, title)
print(p)

# Addendum Figure 2 -- AUC performance of the AI-based healthy index among AI patients and controls
results = readRDS(AI.BASED.HI.IN.PATH)
meta = readRDS(AI.BASED.HI.META.IN.PATH)
title = 'ROC Curve of AI-Based Classifier LOO CV predictions\namong AI and Healthy Subjects'

p = make_auc_plot(results, meta, title)
print(p)

# Addendum Figure 3 -- AUC performance of the original healthy index among PID patients and controls
results = readRDS(ALL.BASED.HI.IN.PATH)
meta = readRDS(ALL.BASED.HI.META.IN.PATH)
results = results[meta$condition %in% c(PID, 'Healthy'), ]
meta = meta[meta$condition %in% c(PID, 'Healthy'), ]
title = 'ROC Curve of Original Classifier LOO CV predictions\namong PID and Healthy Subjects'

p = make_auc_plot(results, meta, title)
print(p)

# Addendum Figure 4 -- AUC performance of the original healthy index among PID patients and controls
results = readRDS(PID.BASED.HI.IN.PATH)
meta = readRDS(PID.BASED.HI.META.IN.PATH)
title = 'ROC Curve of PID-Based Classifier LOO CV predictions\namong PID and Healthy Subjects'

p = make_auc_plot(results, meta, title)
print(p)

# Addendum Figure 5 -- Bar plots of original healthy index for AI subjects
results = readRDS(ALL.BASED.HI.IN.PATH)
meta = readRDS(ALL.BASED.HI.META.IN.PATH)
results = results[meta$condition %in% c(AI, 'Healthy'), ]
meta = meta[meta$condition %in% c(AI, 'Healthy'), ]
title = 'Barplots of Original Classifier LOO CV predictions\namong AI and Healthy Subjects'

p = make_bar_plots(results, meta, title)
print(p)

# Addendum Figure 6 -- Bar plots of AI-based healthy index for AI subjects
results = readRDS(AI.BASED.HI.IN.PATH)
meta = readRDS(AI.BASED.HI.META.IN.PATH)
title = 'Barplots of AI-Based Classifier LOO CV predictions\namong AI and Healthy Subjects'

p = make_bar_plots(results, meta, title)
print(p)

# Addendum Figure 7 -- Bar plots of original healthy index for PID subjects
results = readRDS(ALL.BASED.HI.IN.PATH)
meta = readRDS(ALL.BASED.HI.META.IN.PATH)
results = results[meta$condition %in% c(PID, 'Healthy'), ]
meta = meta[meta$condition %in% c(PID, 'Healthy'), ]
title = 'Barplots of Original Classifier LOO CV predictions\namong PID and Healthy Subjects'

p = make_bar_plots(results, meta, title)
print(p)

# Addendum Figure 8 -- Bar plots of PID-based healthy index for PID subjects
results = readRDS(PID.BASED.HI.IN.PATH)
meta = readRDS(PID.BASED.HI.META.IN.PATH)
title = 'Barplots of PID-Based Classifier LOO CV predictions\namong PID and Healthy Subjects'
p = make_bar_plots(results, meta, title)
print(p)

# Addendum Figure 9 -- p value plots for original classifier
p.vals = readRDS(ALL.BASED.HI.PVALS.IN.PATH)
title = 'Negative log10 adjusted pvalues for feature GVI\nin original classifier'
p = make_pval_plots(p.vals, title = title)
print(p)

# Addendum Figure 10 -- p value plots for the AI-based classifier
p.vals = readRDS(AI.BASED.HI.PVALS.IN.PATH)
title = 'Negative log10 adjusted pvalues for feature GVI\nin AI-based classifier'
p = make_pval_plots(p.vals, title = title)
print(p)

# Addendum Figure 11 -- p value plots for the PID-based classifier
p.vals = readRDS(PID.BASED.HI.PVALS.IN.PATH)
title = 'Negative log10 adjusted pvalues for feature GVI\nin PID-based classifier'
p = make_pval_plots(p.vals, title = title)
print(p)

# Addendum Figure 12 -- correlation between original healthy index and AI-based healthy index
results.all = readRDS(ALL.BASED.HI.IN.PATH)
results.ai = readRDS(AI.BASED.HI.IN.PATH)
meta.all = readRDS(ALL.BASED.HI.META.IN.PATH)
meta.ai = readRDS(AI.BASED.HI.META.IN.PATH)

p = make_scatter_plots_two_group(results.1 = results.all, results.2 = results.ai, 
                                 meta.1 = meta.all, meta.2 = meta.ai, 
                                 title = 'Correlation between original healthy index\nand AI-Based healthy index', 
                                 group1.name = 'Healthy', group2.name = 'AI',
                                 x.label = 'Original HI', y.label = 'AI-Based HI')

print(p)

# Addendum Figure 13 -- correlation between original healthy index and AI-based healthy index
results.all = readRDS(ALL.BASED.HI.IN.PATH)
results.ai = readRDS(AI.BASED.HI.IN.PATH)
meta.all = readRDS(ALL.BASED.HI.META.IN.PATH)
meta.ai = readRDS(AI.BASED.HI.META.IN.PATH)

p = make_scatter_plots_multi_group(results.1 = results.all, results.2 = results.ai, 
                                 meta.1 = meta.all, meta.2 = meta.ai, 
                                 title = 'Correlation between original healthy index\nand AI-Based healthy index',
                                 x.label = 'Original HI', y.label = 'AI-Based HI')
print(p)

# Addendum Figure 14 -- correlation between original healthy index and AI-based healthy index
results.all = readRDS(ALL.BASED.HI.IN.PATH)
results.pid = readRDS(PID.BASED.HI.IN.PATH)
meta.all = readRDS(ALL.BASED.HI.META.IN.PATH)
meta.pid = readRDS(PID.BASED.HI.META.IN.PATH)

p = make_scatter_plots_two_group(results.1 = results.all, results.2 = results.pid, 
                                 meta.1 = meta.all, meta.2 = meta.pid, 
                                 title = 'Correlation between original healthy index\nand PID-Based healthy index', 
                                 group1.name = 'Healthy', group2.name = 'PID',
                                 x.label = 'Original HI', y.label = 'PID-Based HI')
print(p)

# Addendum Figure 15 -- correlation between original healthy index and AI-based healthy index
results.all = readRDS(ALL.BASED.HI.IN.PATH)
results.pid = readRDS(PID.BASED.HI.IN.PATH)
meta.all = readRDS(ALL.BASED.HI.META.IN.PATH)
meta.pid = readRDS(PID.BASED.HI.META.IN.PATH)

p = make_scatter_plots_multi_group(results.1 = results.all, results.2 = results.pid, 
                                 meta.1 = meta.all, meta.2 = meta.pid, 
                                 title = 'Correlation between original healthy index\nand PID-Based healthy index', 
                                 x.label = 'Original HI', y.label = 'PID-Based HI')
print(p)

# Addendum Figure 16 -- correlation between AI-based healthy index and PID-based healthy index among healthy controls
results.ai = readRDS(AI.BASED.HI.IN.PATH)
results.pid = readRDS(PID.BASED.HI.IN.PATH)
meta.ai = readRDS(AI.BASED.HI.META.IN.PATH)
meta.pid = readRDS(PID.BASED.HI.META.IN.PATH)

p = make_scatter_plots_multi_group(results.1 = results.ai, results.2 = results.pid, 
                                   meta.1 = meta.ai, meta.2 = meta.pid, 
                                   title = 'Correlation between AI-Based healthy index\nand PID-Based healthy index among healthy controls', 
                                   x.label = 'AI-Based HI', y.label = 'PID-Based HI',
                                   remove.healthy = F)
print(p)

# Addendum Figures 17 & 18 -- correlation between AI-based healthy index and PID-based healthy index among AI patients
results.pid.from.ai = readRDS(PID.PREDICTIONS.FROM.AI.IN.PATH)
results.pid.from.pid = readRDS(PID.BASED.HI.IN.PATH)
meta.pid = readRDS(PID.BASED.HI.META.IN.PATH)

meta.pid = meta.pid[meta.pid$condition %in% PID,]
p = make_scatter_plots_one_group(results.1 = results.pid.from.ai, results.2 = results.pid.from.pid, 
                                   meta.1 = meta.pid, meta.2 = meta.pid, 
                                   title = 'Correlation between AI-Based healthy index\nand PID-Based healthy index\namong PID subjects', 
                                   x.label = 'AI-Based HI', y.label = 'PID-Based HI')
print(p)

p = make_scatter_plots_multi_group(results.1 = results.pid.from.ai, results.2 = results.pid.from.pid, 
                                 meta.1 = meta.pid, meta.2 = meta.pid, 
                                 title = 'Correlation between AI-Based healthy index\nand PID-Based healthy index\namong PID subjects (condition-specific)', 
                                 x.label = 'AI-Based HI', y.label = 'PID-Based HI')
print(p)

# Addendum Figures 19 & 20 -- correlation between AI-based healthy index and PID-based healthy index among PID patients
results.ai.from.pid = readRDS(AI.PREDICTIONS.FROM.PID.IN.PATH)
results.ai.from.ai = readRDS(AI.BASED.HI.IN.PATH)
meta.ai = readRDS(AI.BASED.HI.META.IN.PATH)

meta.ai = meta.ai[meta.ai$condition %in% AI,]
p = make_scatter_plots_one_group(results.1 = results.ai.from.pid, results.2 = results.ai.from.ai, 
                                 meta.1 = meta.ai, meta.2 = meta.ai, 
                                 title = 'Correlation between PID-Based healthy index\nand AI-Based healthy index\namong AI subjects', 
                                 x.label = 'PID-Based HI', y.label = 'AI-Based HI')
print(p)

p = make_scatter_plots_multi_group(results.1 = results.ai.from.pid, results.2 = results.ai.from.ai, 
                                 meta.1 = meta.ai, meta.2 = meta.ai, 
                                 title = 'Correlation between PID-Based healthy index\nand AI-Based healthy index\namong AI subjects (condition-specific)', 
                                 x.label = 'PID-Based HI', y.label = 'AI-Based HI')
print(p)

# Addendum Figures 21 & 22 -- correlation between AI-based healthy index and original healthy index among PID patients
results.pid.from.ai = readRDS(PID.PREDICTIONS.FROM.AI.IN.PATH)
results.pid.from.all = readRDS(ALL.BASED.HI.IN.PATH)
meta.all = readRDS(ALL.BASED.HI.META.IN.PATH)

meta.all = meta.all[meta.all$condition %in% PID,]
p = make_scatter_plots_one_group(results.1 = results.pid.from.all, results.2 = results.pid.from.ai, 
                                 meta.1 = meta.all, meta.2 = meta.all, 
                                 title = 'Correlation between AI-Based healthy index\nand original healthy index\namong PID subjects', 
                                 x.label = 'Original HI', y.label = 'AI-Based HI')
print(p)

p = make_scatter_plots_multi_group(results.1 = results.pid.from.all, results.2 = results.pid.from.ai, 
                                 meta.1 = meta.all, meta.2 = meta.all, 
                                 title = 'Correlation between AI-Based healthy index\nand original healthy index\namong PID subjects (condition-specific)', 
                                 x.label = 'Original HI', y.label = 'AI-Based HI')
print(p)

# Addendum Figures 23 & 24 -- correlation between AI-based healthy index and original healthy index among AI patients
results.ai.from.pid = readRDS(AI.PREDICTIONS.FROM.PID.IN.PATH)
results.ai.from.all = readRDS(ALL.BASED.HI.IN.PATH)
meta.all = readRDS(ALL.BASED.HI.META.IN.PATH)

meta.all = meta.all[meta.all$condition %in% AI,]
p = make_scatter_plots_one_group(results.1 = results.ai.from.all, results.2 = results.ai.from.pid, 
                                 meta.1 = meta.all, meta.2 = meta.all, 
                                 title = 'Correlation between PID-Based healthy index\nand original healthy index\namong AI subjects', 
                                 x.label = 'Original HI', y.label = 'PID-Based HI')
print(p)

p = make_scatter_plots_multi_group(results.1 = results.ai.from.all, results.2 = results.ai.from.pid, 
                                 meta.1 = meta.all, meta.2 = meta.all, 
                                 title = 'Correlation between PID-Based healthy index\nand original healthy index\namong AI subjects (condition-specific)', 
                                 x.label = 'Original HI', y.label = 'PID-Based HI')
print(p)
dev.off()
