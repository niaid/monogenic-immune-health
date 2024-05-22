library(tidyverse)
library(cowplot)
library(ggraph)
library(tidygraph)
library(ggpubr)
library(ggrepel)

JIVE.PC.IN.PATH <- snakemake@input[["jive_pcs"]]#"Integration_output/jive/subject/prcomp_list.rds"
JIVE.IN.PATH <- snakemake@input[["jive"]]#"Integration_output/jive/subject/jive.rds"

BOXPLOT.OUT.PATH <- snakemake@output[["boxplots"]]#"Paper_1_Figures/Supplemental_Figure_3/jive_array_indiv_boxplots.pdf"
SCATTER.OUT.PATH <- snakemake@output[["scatter"]]#"Paper_1_Figures/Supplemental_Figure_3/jive_array_indiv_scatter_w_centroids.pdf"
source("scripts/util/paper/abbrev_cond.R")

prcomp.list <- readRDS(JIVE.PC.IN.PATH)
jive <- readRDS(JIVE.IN.PATH)
pdat <- jive$pdat
array_indiv <- prcomp.list$array.ind$x

stopifnot(identical(rownames(array_indiv), pdat$patient_id))

array_indiv <- array_indiv %>% 
  as.data.frame() %>%
  bind_cols(pdat) %>% 
  mutate(cond.abbrev = abbrev_cond(condition)) %>%
  mutate(cond.grouped = group_cond(condition))


#Inspired by this post
#https://stackoverflow.com/questions/23463324/r-add-centroids-to-scatter-plot

array_indiv_centroids <- array_indiv %>% group_by(condition) %>%
        summarise(mean_pc1 = mean(PC1), mean_pc2 = mean(PC2),n_subj = n())

array_indiv <- left_join(array_indiv, array_indiv_centroids)

array_indiv_sub <- array_indiv %>% 
        filter(n_subj > 3) %>%
        select(mean_pc1, mean_pc2, cond.abbrev, cond.grouped) %>% distinct()

pca.plot.points <- 
  ggplot(array_indiv, aes(x = PC1, y = PC2, color = cond.grouped)) +
  #geom_point( size = 3) +
  geom_text(aes(label = cond.abbrev), size = 2) +
  geom_point(data = array_indiv_sub, aes(x=mean_pc1, y=mean_pc2),size=5)+
  geom_segment(aes(x=mean_pc1, y=mean_pc2, xend=PC1, yend=PC2), alpha = .2) +
  #geom_text_repel(data = array_indiv_sub, aes(x=mean_pc1, y=mean_pc2, label = cond.abbrev),size=5, color = "black")+
  geom_text_repel(data = array_indiv_sub, aes(x=mean_pc1, y=mean_pc2, label = cond.abbrev),size=5)+
  theme_bw() +
  theme(legend.position = "none") #+
  #theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  #theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

pdf(SCATTER.OUT.PATH, height = 5, width = 6)
print(pca.plot.points)
dev.off()

pc.medians <-
  array_indiv %>%
  group_by(cond.abbrev) %>%
  summarise(pc1.median = median(PC1), pc2.median = median(PC2))


pc1.order <- pc.medians$cond.abbrev[order(pc.medians$pc1.median)]
pc1.order <- c("Healthy", pc1.order[pc1.order != "Healthy"])
array_indiv$cond.abbrev <- factor(array_indiv$cond.abbrev, levels = pc1.order)

pc1.box <- 
  ggplot(array_indiv, aes(x = cond.abbrev, y = PC1)) +
  geom_boxplot(outlier.shape = NA, aes(fill = cond.grouped)) +
  ggbeeswarm::geom_beeswarm(size = .8, alpha = .4)+
  theme_bw() +
  stat_compare_means(ref.group = "Healthy", hide.ns = TRUE, label = "p.signif", color = "red") +
  coord_flip() + 
  geom_vline(xintercept = 1.5) +
  theme(legend.position = "none") +
  xlab("Condition") +
  ylab("iPC1")


pc2.order <- pc.medians$cond.abbrev[order(pc.medians$pc2.median)]
pc2.order <- c("Healthy", pc2.order[pc2.order != "Healthy"])
array_indiv$cond.abbrev <- factor(array_indiv$cond.abbrev, levels = pc2.order)

pc2.box <- 
  ggplot(array_indiv, aes(x = cond.abbrev, y = PC2)) +
  geom_boxplot(outlier.shape = NA, aes(fill = cond.grouped)) +
  ggbeeswarm::geom_beeswarm(size = .8, alpha = .4)+
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip() + 
  stat_compare_means(ref.group = "Healthy", hide.ns = TRUE, label = "p.signif", color = "red") +
  theme(legend.position = "none") +
  geom_vline(xintercept = 1.5) +
  xlab("Condition") +
  ylab("iPC2")

pdf(BOXPLOT.OUT.PATH, height = 3, width = 4)
print(pc1.box)
print(pc2.box)
dev.off()
