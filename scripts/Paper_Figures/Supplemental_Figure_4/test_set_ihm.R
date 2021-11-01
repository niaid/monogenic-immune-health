library(tidyverse)
library(ggpubr)

setwd("../../..")

source('scripts/util/Plotting/plot_auc.R')

PRED.IN.PATH <- "Pipeline_out/Classification/predictions/healthy_rf_testing_predictions_all.RDS"

META.IN.PATH <- "Pipeline_out/Classification/meta_data/healthy_random_forest_testing_sample_meta_data_all.RDS"

ROC.PLOT.OUT.PATH <- "Pipeline_out/Paper_1_Figures_20210929/Supplemental_Figure_4/test_set_predictions_roc.pdf"

BOX.PLOT.OUT.PATH <- "Pipeline_out/Paper_1_Figures_20210929/Supplemental_Figure_4/test_set_predictions_boxplot.pdf"

results <- readRDS(PRED.IN.PATH)

meta <- readRDS(META.IN.PATH)

stopifnot(identical(rownames(results), meta$patient_id))

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
ggsave(ROC.PLOT.OUT.PATH, p, device = 'pdf', height = 2, width = 2)

plot_dat <- meta %>%
        mutate(IHM = result) %>%
        mutate(cond2 = ifelse(condition == "Healthy", "Healthy", "Disease"))

p <- ggplot(plot_dat, aes(x = cond2, y = IHM)) +
        geom_boxplot(outlier.shape = NA) +
        #geom_text(aes(label = condition, color = condition)) +
        ggbeeswarm::geom_beeswarm(aes(color = condition)) +
        xlab("") +
        stat_compare_means()

ggsave(BOX.PLOT.OUT.PATH, p, device = 'pdf', height = 4, width = 4)
