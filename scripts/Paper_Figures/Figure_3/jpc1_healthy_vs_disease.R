library(ggplot2)
library(ggpubr)
library(dplyr)

JIVE.PC.IN.PATH <- snakemake@input[["jive_pcs"]]#"Integration_output/jive/subject/prcomp_list.rds"
JIVE.IN.PATH <- snakemake@input[["jive"]]#"Integration_output/jive/subject/jive.rds"

FIG.OUT.PATH <- snakemake@output[[1]]#"Paper_1_Figures/Figure_3/jpc1_hc_vs_all_try_colorschemes.pdf"
#Wilcoxon of Healthy vs Not Healthy pc1

source("scripts/util/paper/abbrev_cond.R")

prcomp.list <- readRDS(JIVE.PC.IN.PATH)
jive <- readRDS(JIVE.IN.PATH)
pdat <- jive$pdat
joint <- prcomp.list$joint$x

stopifnot(identical(rownames(joint), pdat$patient_id))

joint <- joint %>% 
  as.data.frame() %>%
  bind_cols(pdat) %>% 
  mutate(cond.abbrev = abbrev_cond(condition)) %>%
  mutate(cond.grouped = group_cond(condition))

#joint <- joint %>% mutate(Healthy = condition == "Healthy")
joint <- joint %>% mutate(Healthy = ifelse(condition == "Healthy", "Healthy", "Disease"))
#file.remove("Paper_1_Figures/Figure_3/jpc1_hc_vs_all.pdf")

p <- ggplot(joint, aes(x = Healthy, y = PC1)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_bw() +
  stat_compare_means() + xlab("") + ylab("jPC1")
ggsave(plot = p, filename = FIG.OUT.PATH, height = 3, width = 3)
