
library(tidyverse)
library(ggpubr)

IHM.IN.PATH <- snakemake@input[[1]]
META.IN.PATH <- snakemake@input[[2]]

FIG.OUT.PATH <- snakemake@output[[1]]

ihm <- readRDS(IHM.IN.PATH)

meta <- readRDS(META.IN.PATH)

#mods <- readRDS("../Pipeline_out/Classification/results_no_pm2/healthy_rf_models_all.RDS")
#
#hmod <- mods$all.modules.plus.grey.with.tbnks
#
##don't see soma purple mod
#hmod$importance

ihm <- ihm$all.modules.plus.grey.with.tbnks

dat <- meta %>%
        mutate(ihm = ihm)


#keep_cond <- c("47CGD", "XCGD", "Healthy", "Job", "STAT1 GOF", "FMF")
keep_cond <- c("Healthy")
p <- dat %>%
        filter(condition %in% keep_cond) %>%
        ggplot(aes(x = Age, y = ihm)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE) +
        stat_cor(method = "spearman") +
        theme_bw() +
        facet_wrap(~condition)

ggsave(plot = p, filename = FIG.OUT.PATH, height = 3, width = 3)

