library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggbeeswarm)

if(!exists("snakemake")){
  setwd("../../..")
  source("scripts/util/paper/parse_snakemake.R")
  parse_snakemake(rule = "supplemental_figure_3_cgd_jpc1")
}

JIVE.PC.IN.PATH <- snakemake@input[["jive_pcs"]]
JIVE.IN.PATH <- snakemake@input[["jive"]]

FIG.OUT.PATH <- snakemake@output[[1]]

prcomp.list <- readRDS(JIVE.PC.IN.PATH)
jive <- readRDS(JIVE.IN.PATH)
pdat <- jive$pdat
joint <- prcomp.list$joint$x

stopifnot(identical(rownames(joint), pdat$patient_id))

joint <- joint %>% 
  as.data.frame() %>%
  bind_cols(pdat)

joint <- joint %>% filter(grepl("CGD", condition))

pdf(FIG.OUT.PATH, height = 3, width =3)
p <- ggplot(joint, aes(x = condition, y = PC1)) +
        geom_boxplot(outlier.shape = NA) +
        geom_beeswarm() +
        stat_compare_means() +
        ylab("jPC1") + xlab("")
print(p)


p <- ggplot(joint, aes(x = condition, y = PC1)) +
        geom_boxplot(outlier.shape = NA) +
        geom_beeswarm() +
        stat_compare_means(method = "t.test") +
        ylab("jPC1") + xlab("")
print(p)
dev.off()

#ggsave(plot = p, filename = FIG.OUT.PATH, height = 3, width = 3)
