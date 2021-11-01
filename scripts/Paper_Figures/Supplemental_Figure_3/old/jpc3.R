suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(ggraph)
  library(tidygraph)
})

prcomp.list <- readRDS("Integration_output/jive/subject/prcomp_list.rds")
source("scripts/util/JIVE/JIVE_extract.R")
source("scripts/util/Trees/mst_related.R")
source("scripts/util/paper/abbrev_cond.R")
source("scripts/util/Permutation/perm_2group.R")

joint <- prcomp.list$joint$x
stopifnot(identical(rownames(joint), prcomp.list$pdat$patient_id))

joint <- joint %>% as.data.frame() %>%
  bind_cols(prcomp.list$pdat) %>% 
  mutate(cond.abbrev = abbrev_cond(condition)) %>%
  mutate(cond.grouped = group_cond(condition))

pc.medians <-
  joint %>%
  group_by(cond.abbrev) %>%
  summarise(pc3.median = median(PC3))


pc3.order <- pc.medians$cond.abbrev[order(-pc.medians$pc3.median)]
joint$cond.abbrev <- factor(joint$cond.abbrev, levels = pc3.order)

#keep same color scheme as other plots
pc3.box <- 
  ggplot(joint, aes(x = cond.abbrev, y = PC3, color = cond.grouped)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_bw() +
  coord_flip() + 
  theme(axis.title.x.bottom = element_blank(), 
        axis.title.y = element_blank(), 
        legend.position = "none")

pc3.box
