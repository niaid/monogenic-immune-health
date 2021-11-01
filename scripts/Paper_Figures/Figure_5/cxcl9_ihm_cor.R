library(tidyverse)
library(Biobase)
library(ggpubr)

HI.IN.PATH = snakemake@input[["rf_results"]]#'Classification/results/healthy_rf_results_all.RDS'
META.IN.PATH = snakemake@input[["rf_meta"]]#'Classification/random_forest_sample_meta_data.RDS'
SOMA.IN.PATH = snakemake@input[["soma_data"]]#'Data/Somalogic/data_analysis_ready/analysis_ready_subject_level_training_somalogic.rds'

FIG.OUT.PATH = snakemake@output[[1]]#"Paper_1_Figures/Figure_5/cxcl9_cor_immune_health_metric.pdf"

hi_dat <- readRDS(HI.IN.PATH)
meta <- readRDS(META.IN.PATH)
soma <- readRDS(SOMA.IN.PATH)


hi_dat <- hi_dat %>%
        rownames_to_column(var = "patient_id") %>%
        mutate(`Immune Health Metric` = all.modules.plus.grey.with.tbnks) %>%
        select(patient_id, `Immune Health Metric`)

soma_mat <- exprs(soma)
grep("9", rownames(soma_mat), ignore.case = T, value = T)
grep("mig", rownames(soma_mat), ignore.case = T, value = T)

feat_dat <- featureData(soma)@data

feat_dat %>% filter(Target == "MIG")

cxcl9_dat <- pData(soma) %>%
        select(patient_id, condition) %>%
        mutate(`CXCL9` = soma_mat["MIG", ])

plot_dat <- left_join(hi_dat, cxcl9_dat)

plot_dat <- plot_dat %>%
        mutate(healthy = ifelse(condition == "Healthy", "Healthy", "Disease")) %>%
        mutate(healthy = factor(healthy, levels = c("Healthy", "Disease")))

p <- ggplot(plot_dat, aes(x = `CXCL9`, y = `Immune Health Metric`)) +
        geom_point() +
        stat_cor(method = "spearman") +
        theme_bw() +
        facet_wrap(~healthy) +
        ggtitle("Monogenic data using Immune Health Metric scores directly")

ggsave(plot = p, filename = FIG.OUT.PATH, height = 3, width = 6)

#p <- ggplot(plot_dat %>% filter(condition != "Healthy"), aes(x = `CXCL9`, y = `Immune Health Metric`)) +
#        geom_text(aes(label = condition, color= condition)) +
#        stat_cor(method = "spearman") +
#        theme_bw() +
#        #facet_wrap(~healthy) +
#        ggtitle("Monogenic data using Immune Health Metric scores directly")
#
#ggsave(plot = p, filename = "Misc/cxcl9_cor_immune_health_metric_healthy_contols_separate_color_condition.pdf", height = 8, width = 15)

