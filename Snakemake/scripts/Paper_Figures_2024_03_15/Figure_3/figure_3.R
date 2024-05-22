# Load libraries
library(ggplot2)
library(reshape2)
library(cowplot)
library(ggpubr)
library(dplyr)
library(tidyr)
library(ggpubr)

# Set paths
HI.IN.PATH = snakemake@input[[1]]#'Classification/results/healthy_rf_results_all.RDS'
META.IN.PATH = snakemake@input[[2]]#'Classification/random_forest_sample_meta_data.RDS'
PC1.IN.PATH = snakemake@input[[3]]#"Integration_output/jive/subject/prcomp_list.rds"
FEATURE.GVI.PVALS.IN.PATH = snakemake@input[[4]]#"Classification/results/healthy_rf_pvals_all.RDS"
TRANSCRIPTIONAL.SURROGATE.SIGNATURE.ENRICHMENTS.IN.PATH = snakemake@input[[5]]#"Classification/transcriptional_surrogates/surrogate_enrichments.RDS"
DESIGN.MAT.IN.PATH <- snakemake@input[[6]]

FIGURE.4b.OUT.PATH = snakemake@output[[1]]#"Paper_1_Figures/Figure_4/4.b.pdf",
FIGURE.4cde.OUT.PATH = snakemake@output[[2]]#"Paper_1_Figures/Figure_4/4.cde.pdf"
FIGURE.4f.OUT.PATH = snakemake@output[[3]]#"Paper_1_Figures/Figure_4/4.f.pdf"
FIGURE.4g.OUT.PATH = snakemake@output[[4]]#"Paper_1_Figures/Figure_4/4.g.pdf"
FIGURE.4h.OUT.PATH = snakemake@output[[5]]#"Paper_1_Figures/Figure_4/4.h.pdf"

# We source utilities
source('scripts/util/Plotting/enrichments.R')
source('scripts/util/Plotting/plot_auc.R')
source('scripts/util/paper/abbrev_cond.R')
source('scripts/util/Groups/groups.R')

# Figure 4a -- Cartoon displaying classifier schema

# Figure 4b -- AUC curve for classifier

## Load the healthy rf prediction results and the sample metadata
results = readRDS(HI.IN.PATH)
meta = readRDS(META.IN.PATH)

## Get the "all-features" classifier predictions, which we have been using for the Healthy Index
result = results$all.modules.plus.grey.with.tbnks

## Create the ROC curve from the HI
roc = get_roc(result, meta$condition, 'Healthy')

## Get the AUC associated with this ROC curve
auc = get_auc(roc)

## Round the AUC to 2 digits
auc = format(auc, digits = 2)

## Plot the ROC curve with the AUC displayed
p = ggplot(roc, aes(x = fpr, y = tpr)) + geom_line(color = 'black', show.legend = FALSE) + 
  geom_abline(slope = 1, linetype = 'dashed', color = 'grey') + 
  theme_bw() + geom_text(aes(x = .75, y = .25), size = 4, label = paste0('AUC: ', auc), show.legend = FALSE) + 
  xlab('False Positive Rate') + ylab('True Positive Rate') + 
  theme(axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8)
  )

## Save plot
ggsave(FIGURE.4b.OUT.PATH, p, device = 'pdf', height = 2, width = 2)

# Figure 4c -- healthy index barplots for each condition, arranged by condition supertype

## Load the healthy index and subject meta-data
results = readRDS(HI.IN.PATH)
meta = readRDS(META.IN.PATH)

## Get the conditions associated with each supertype
AI = util.get_ai()
PID = c(util.get_pid(), 'NEMO') # We put in NEMO manually because it shows up as non-PID in the database
Telo = util.get_tert_terc()

## Create a data frame with the HI and condition for each subject
df = data.frame(healthy.index = results$all.modules.plus.grey.with.tbnks, 
                condition = meta$condition) %>%
  mutate(condition = as.character(condition)) %>%
  mutate(group = condition %>% # Add in condition super-type
           replace(condition %in% AI, 'AI') %>%
           replace(condition %in% PID, 'PID') %>%
           replace(condition %in% Telo, 'Telo')) %>%
  mutate(group = factor(group, levels = c('Healthy','AI','Telo','PID'))) %>%
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

## Plot the box plots
HI_max = max(df$healthy.index) + .01
p1 = ggplot(df, aes(x = condition, y = healthy.index, fill = group)) + 
  geom_boxplot(outlier.colour = NA) + 
  ylim(0, HI_max) +
  theme_bw() + 
  geom_jitter() + 
  coord_flip() + 
  stat_compare_means(ref.group = "Healthy", hide.ns = TRUE, label = "p.signif", color = "red", label.y = 0, size = 10) +
  theme(axis.text.x = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.key.size = unit(2,"line"))

# Figure 4d -- density plot for major condition groups

## Load data
results = readRDS(HI.IN.PATH)
meta = readRDS(META.IN.PATH)

## Choose the conditions to plot
conditions = c('Healthy','47CGD','XCGD','Job','STAT1 GOF','FMF')

## Create a data frame with the HI and condition for each subject from the conditions of interest
df = data.frame(healthy.index = results$all.modules.plus.grey.with.tbnks, 
                condition = meta$condition) %>%
  mutate(condition = as.character(condition)) %>%
  filter(condition %in% conditions) %>%
  mutate(condition = abbrev_cond(condition)) %>%
  mutate(condition = factor(condition, levels = abbrev_cond(conditions)))

## Make the density plots
p2 = ggplot(df, aes(x = healthy.index, fill = condition)) + geom_density(alpha = .4) +
  ylab('Density') + xlab('Healthy Index') + xlim(0, HI_max) + 
  theme_bw() + theme(
    axis.text.x = element_text(size = 15),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15)
  )

# Figure 4e -- correlation between healthy index and PC1 score in healthy and conditions

## Load Joint PC1
jive = readRDS(PC1.IN.PATH)
PC1 = jive$joint$x[,1]

## Load healthy index
results = readRDS(HI.IN.PATH)
predictions = results$all.modules.plus.grey.with.tbnks
names(predictions) = rownames(results)

## Load patient metadata
meta = readRDS(META.IN.PATH)

## Get the patients that have a healthy index and joint PC1
ids = intersect(names(predictions), names(PC1))
PC1 = PC1[ids]
predictions = predictions[ids]
conditions = meta[ids, 'condition']
group = ifelse(conditions == "Healthy", "Healthy", "Disease")
group = factor(group, levels = c("Healthy", "Disease"))

## Create a data frame joining the two
df = data.frame(predictions, PC1, group)

## Create a scatter plot with regression lines showing the relationships between HI and PC1
p3 = ggplot(df, aes(x = predictions, y = PC1, color = group)) + ylab('PC1 Score') +
  xlab('Healthy Index') +
  geom_point() + geom_smooth(method = 'lm', formula = y ~ x, se = FALSE) +
  ylim(-50, 50) + xlim(0, HI_max) +
  stat_cor(label.y = c(30,45), label.x = c(0,0), show.legend = FALSE, size = 5) +
  theme_bw() + 
  theme(axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))

## Plot Figure 4c-e
p = plot_grid(plotlist = list(p1, p2, p3), nrow = 3, align = 'v', rel_heights = c(10, 3, 6))
ggsave(FIGURE.4cde.OUT.PATH, p, height = 10, width = 10)

# Figure 4f -- Age relationships of top conditions 

## Load data
results = readRDS(HI.IN.PATH)
meta = readRDS(META.IN.PATH)

## Choose the conditions to plot
conditions = c('Healthy','47CGD','XCGD','Job','STAT1 GOF','FMF')

## Make a data frame with the healthy indexes, conditions, and ages of the patients from the conditions of interest

library(tidyverse)
df <- results %>%
        rownames_to_column(var = "patient_id") %>%
        select(patient_id, all.modules.plus.grey.with.tbnks) %>%
        rename(healthy.index = all.modules.plus.grey.with.tbnks) %>%
        left_join(meta) %>%
  filter(condition %in% conditions) %>%
  mutate(condition = as.character(condition)) %>%
  mutate(condition = abbrev_cond(condition)) %>%
  mutate(condition = factor(condition, levels = abbrev_cond(conditions)))

p = ggplot(df, aes(x = Age, y = healthy.index)) + geom_point() + 
  facet_wrap(~condition, ncol = 3, nrow = 2) + ylab('Healthy Index') +
  geom_smooth(method = 'lm', formula = y ~ x, se = FALSE) +
  ylim(0, .6) + xlim(0, 80) +
  stat_cor(label.x = 0, label.y = .55) +
  theme_bw() + 
  theme(axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        strip.text = element_text(size = 15))

ggsave(FIGURE.4f.OUT.PATH, p, height = 5, width = 10)

# Figure 4g -- Bar charts of top features of the classifier

## Load the HI feature pvalues
results = readRDS(FEATURE.GVI.PVALS.IN.PATH)
p.vals = results$all.modules.plus.grey.with.tbnks
meta = readRDS(META.IN.PATH)

design_mat_list <- readRDS(DESIGN.MAT.IN.PATH)
design_mat <- design_mat_list$all.modules.plus.grey.with.tbnks

## Adjust the p-values using an FDR and convert to negative log10 pvalues
p.adjusted = p.adjust(p.vals, 'fdr')
neg.log10.p.adjusted = -1 * log10(p.adjusted)

## Create a data frame with the feature names, data type, and negative log 10 pvalues, just for features passing the FDR cutoff
df = data.frame(label = names(p.vals), p.adjusted = p.adjusted, neg.log10.pvals = neg.log10.p.adjusted) %>%
  filter(p.adjusted < .20) %>%
  select(-p.adjusted)

keep_feat <- df$label %>% as.character()

feat_mat <- design_mat[, keep_feat] 

feat_dat <- feat_mat %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "patient_id") %>%
        gather(key = feat, value = value, -patient_id) %>%
        left_join(meta) %>%
        mutate(healthy = condition == "Healthy")

ttestf <- function(x) {
        t.test(value ~ healthy, paired = FALSE, data = x) %>% broom::tidy()
}
t_test_dat <- feat_dat %>%
        group_by(feat) %>%
        do(ttestf(.)) %>% ungroup()

feat_sign <- t_test_dat %>%
        mutate(t_sign = sign(statistic)) %>%
        select(feat, t_sign) %>%
        tibble::deframe()



df <- df %>%
        mutate(t_sign = feat_sign[match(label, names(feat_sign))])

df = df %>%
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
           gsub('nk_cells_percent','NK Cells (%)', .) %>%
           gsub('nk_cells_abs','NK Cells (#)', .) %>%
           gsub('MIP.1a','MIP 1a', .) %>%
           gsub('purple','PM2', .) %>%
           gsub('Cathepsin.H','Cathepsin H', .) %>%
           gsub('IL.18.Ra','IL-18 Ra', .) %>%
           gsub('rdw','RDW', .) %>%
           gsub('LD78.beta','LD78 b', .))

## We order the features by negative log 10 pvalue
df = df %>% 
  arrange(neg.log10.pvals) %>%
  mutate(label = factor(label, levels = label)) %>%
  mutate(data.type = factor(data.type))

## We plot the feature p-values in a bar plot
p = ggplot(df, aes(y = neg.log10.pvals * t_sign, x = label, fill = data.type)) + 
  geom_bar(stat="identity") + theme_bw() + xlab('Parameter') + ylab('Negative log10 q-values') + coord_flip() +
  #scale_fill_manual(values = c('darkblue','steelblue', 'lightblue')) + 
  labs(fill = 'Data Type') +
  geom_hline(aes(yintercept = -log10(.20), linetype = 'FDR = .20'), color = 'gray', size = 1) +
  geom_hline(aes(yintercept = log10(.20), linetype = 'FDR = .20'), color = 'gray', size = 1) +
  scale_linetype_manual(values = 'dashed') +
  scale_fill_manual(values = c("seagreen4", "violetred", "plum1")) +
  facet_grid(t_sign~1, scales = "free_y", space = "free_y") +
  ylab("-log10(q) * direction") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        strip.background = element_blank(),
        strip.text = element_blank()
  )

## We save the plot
ggsave(FIGURE.4g.OUT.PATH, p, height = 3, width = 5)

# Figure 4h -- barplots of enrichments for gene surrogate signatures of proteins
## We load the enrichments on the surrogate signatures
enrichments = readRDS(TRANSCRIPTIONAL.SURROGATE.SIGNATURE.ENRICHMENTS.IN.PATH)

## We use the enrichment bar plot utility function to create a barplot for two proteins of interest: SAA and MIP1a
p1 = make_enrichment_bar_plot(enrichment = enrichments$somalogic.modules.purple$positive) + 
  ggtitle('Proteomic Purple Module') + scale_fill_manual(values = c('tomato2', 'forestgreen', 'darkorchid2')) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        plot.title = element_text(size = 15))

p2 = make_enrichment_bar_plot(enrichment = enrichments$somalogic.grey.SAA$positive) +
  ggtitle('SAA') +
  scale_fill_manual(values = c('tomato2', 'darkorchid2')) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        plot.title = element_text(size = 15))

## We put the plots together
p = plot_grid(plotlist = list(p1, p2), nrow = 2, align = 'v')

## And save them
ggsave(FIGURE.4h.OUT.PATH, p, height = 6, width = 10)
