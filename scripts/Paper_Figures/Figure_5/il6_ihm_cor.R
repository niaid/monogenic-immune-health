library(tidyverse)
library(Biobase)
library(ggpubr)

#if(exists("snakemake")){
HI.IN.PATH = snakemake@input[["rf_results"]]#'Classification/results/healthy_rf_results_all.RDS'
META.IN.PATH = snakemake@input[["rf_meta"]]#'Classification/random_forest_sample_meta_data.RDS'
SOMA.IN.PATH = snakemake@input[["soma_data"]]#'Data/Somalogic/data_analysis_ready/analysis_ready_subject_level_training_somalogic.rds'

PLOT1.OUT.PATH = snakemake@output[[1]]#"Paper_1_Figures/Figure_5/il6_cor_immune_health_metric.pdf"
PLOT2.OUT.PATH = snakemake@output[[2]]#"Paper_1_Figures/Figure_5/il6_cor_immune_health_metric.pdf"
#} else {
#  HI.IN.PATH = 'Classification/results/healthy_rf_results_all.RDS'
#  META.IN.PATH = 'Classification/random_forest_sample_meta_data.RDS'
#  SOMA.IN.PATH = 'Data/Somalogic/data_analysis_ready/analysis_ready_subject_level_training_somalogic.rds'
#
#  PLOT.OUT.PATH = "Paper_1_Figures/Figure_5/il6_cor_immune_health_metric.pdf"
#
#  setwd("/hpcdata/sg/sg_data/PROJECTS/Monogenic_Project")
#}

source("scripts/util/paper/abbrev_cond.R")

hi_dat <- readRDS(HI.IN.PATH)
meta <- readRDS(META.IN.PATH)
soma <- readRDS(SOMA.IN.PATH)


hi_dat <- hi_dat %>%
        rownames_to_column(var = "patient_id") %>%
        mutate(`Immune Health Metric` = all.modules.plus.grey.with.tbnks) %>%
        select(patient_id, `Immune Health Metric`)

soma_mat <- exprs(soma)
#grep("il.6", rownames(soma_mat), ignore.case = T, value = T)

il6_dat <- pData(soma) %>%
        select(patient_id, condition) %>%
        mutate(`IL-6` = soma_mat["IL.6", ])

plot_dat <- left_join(hi_dat, il6_dat)

plot_dat <- plot_dat %>% 
        mutate(condition = abbrev_cond(condition))

plot_dat <- plot_dat %>%
        mutate(healthy = ifelse(condition == "Healthy", "Healthy", "Disease")) %>%
        mutate(healthy = factor(healthy, levels = c("Healthy", "Disease")))


p <- ggplot(plot_dat, aes(x = `IL-6`, y = `Immune Health Metric`)) +
        geom_point() +
        stat_cor(method = "spearman") +
        theme_bw() +
        facet_wrap(~healthy) +
        ggtitle("Monogenic data using Immune Health Metric scores directly")

ggsave(plot = p, filename = PLOT1.OUT.PATH, height = 3, width = 6)


p <- ggplot(plot_dat %>% filter(condition != "Healthy"), aes(x = `IL-6`, y = `Immune Health Metric`)) +
        geom_text(aes(label = condition, color= condition)) +
        stat_cor(method = "spearman") +
        theme_bw() +
        #facet_wrap(~healthy) +
        ggtitle("Monogenic data using Immune Health Metric scores directly")

ggsave(plot = p, filename = PLOT2.OUT.PATH, height = 8, width = 15)

