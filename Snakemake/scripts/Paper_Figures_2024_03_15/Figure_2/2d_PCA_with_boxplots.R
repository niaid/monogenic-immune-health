library(tidyverse)
library(cowplot)
library(ggraph)
library(tidygraph)
library(ggpubr)
library(ggrepel)

JIVE.PC.IN.PATH <- snakemake@input[["jive_pcs"]]#"Integration_output/jive/subject/prcomp_list.rds"
JIVE.IN.PATH <- snakemake@input[["jive"]]#"Integration_output/jive/subject/jive.rds"
COND.GROUPS.IN.PATH <- snakemake@input[["cond_groups"]]#"Reference/condition_groups.rds"

JPC.FIG.OUT.PATH <- snakemake@output[["jpc_scatter_box"]]#"Paper_1_Figures/Figure_3/jive_scatter_plus_boxplots_colorschemes.pdf"
MAD.PLOT.OUT.PATH <- snakemake@output[["mad_plot"]]#"Paper_1_Figures/Figure_3/mad_plot.pdf"

source("scripts/util/paper/abbrev_cond.R")

prcomp.list <- readRDS(JIVE.PC.IN.PATH)
jive <- readRDS(JIVE.IN.PATH)
pdat <- jive$pdat
joint <- prcomp.list$joint$x

disease_cat <- readRDS(COND.GROUPS.IN.PATH)

stopifnot(identical(rownames(joint), pdat$patient_id))

joint <- joint %>% 
  as.data.frame() %>%
  bind_cols(pdat) %>% 
  mutate(cond.abbrev = abbrev_cond(condition)) %>%
  mutate(cond.grouped = group_cond(condition))


#Inspired by this post
#https://stackoverflow.com/questions/23463324/r-add-centroids-to-scatter-plot

joint_centroids <- joint %>% group_by(condition) %>%
        summarise(mean_pc1 = mean(PC1), mean_pc2 = mean(PC2),
                  med_pc1 = median(PC1), med_pc2 = median(PC2),
                  mad_pc1 = mad(PC1), mad_pc2 = mad(PC2),
                  sd_pc1 = sd(PC1), sd_pc2 = sd(PC2),
                  sem_pc1 = sd(PC1) / sqrt(n()), sem_pc2 = sd(PC2) / sqrt(n()),
                  n = n())

joint <- left_join(joint, joint_centroids)

joint <- left_join(joint, disease_cat)
joint <- joint %>%
        mutate(cond_group = as.character(cond_group)) %>%
        mutate(cond_group = factor(cond_group, 
                                   levels = c("Healthy", "AI", "TERT.TERC", "PID"))) %>%
        mutate(cond_group = replace(cond_group, condition == "Healthy", "Healthy"))

joint_sub <- joint %>% select(mean_pc1, mean_pc2, cond.abbrev, sd_pc1, sd_pc2, cond_group, n, mad_pc1, mad_pc2) %>% distinct()

joint_sub <- joint_sub %>%
        filter(n > 3)

pca.plot.points.connected.ellipse <- 
  ggplot(data = joint, aes(color = cond_group, fill= cond_group)) +
  geom_text(aes(x = PC1, y = PC2, label = cond.abbrev), size = 2.5) +
  geom_point(data = joint_sub, aes(x=mean_pc1, y=mean_pc2),size=5)+
  geom_segment(aes(x=mean_pc1, y=mean_pc2, xend=PC1, yend=PC2), alpha = .2) +
  geom_text_repel(data = joint_sub, aes(x=mean_pc1, y=mean_pc2, label = cond.abbrev),size=5)+
  theme_bw() +
  xlim(c(-47, 40)) +
  ylim(-25, 50) +
  theme(legend.position = "none") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

pc.medians <-
  joint %>%
  group_by(cond.abbrev) %>%
  summarise(pc1.median = median(PC1), pc2.median = median(PC2))




#pc1.order <- pc.medians$cond.abbrev[order(-pc.medians$pc1.median)]
pc1.order <- pc.medians$cond.abbrev[order(pc.medians$pc1.median)]
pc1.order <- c("Healthy", pc1.order[pc1.order != "Healthy"])
joint$cond.abbrev <- factor(joint$cond.abbrev, levels = pc1.order)

#keep same color scheme as other plots
#pc1.box <- 
#  ggplot(joint, aes(x = cond.abbrev, y = PC1, color = cond.grouped)) +
#  geom_boxplot(outlier.shape = NA) +
#  geom_jitter() + 
#  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#  theme_bw() +
#  coord_flip() + 
#  theme(axis.title.x.bottom = element_blank(), 
#        axis.title.y = element_blank(), 
#        legend.position = "none")

pc1.box <- 
  ggplot(joint, aes(x = cond.abbrev, y = PC1)) +
  geom_boxplot(outlier.shape = NA, aes(fill = cond_group)) +
  geom_jitter() + 
  theme_bw() +
  coord_flip() + 
  ylim(c(-47, 40)) +
  theme(legend.position = "none") +
  stat_compare_means(ref.group = "Healthy", hide.ns = TRUE, label = "p.signif", color = "red", label.y = -45) +
  geom_vline(xintercept = 1.5) +
  xlab("Condition")

#pc1.box

pc2.order <- pc.medians$cond.abbrev[order(pc.medians$pc2.median)]
pc2.order <- c("Healthy", pc2.order[pc2.order != "Healthy"])
joint$cond.abbrev <- factor(joint$cond.abbrev, levels = pc2.order)

#pc2.box <- 
#  ggplot(joint, aes(x = cond.abbrev, y = PC2, color = cond.grouped)) +
#  geom_boxplot(outlier.shape = NA) +
#  geom_jitter() + 
#  theme_bw() +
#  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#  #scale_x_discrete(position = "top") +
#  theme(axis.title.y.left = element_blank(), 
#        axis.title.x = element_blank(), legend.position = "none")

pc2.box <- 
  ggplot(joint, aes(x = cond.abbrev, y = PC2)) +
  geom_boxplot(outlier.shape = NA, aes(fill = cond_group)) +
  geom_jitter() + 
  ylim(-25, 50) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .4)) +
  #coord_flip() + 
  theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank()) +
  stat_compare_means(ref.group = "Healthy", hide.ns = TRUE, label = "p.signif", color = "red", label.y = -23) +
  geom_vline(xintercept = 1.5) +
  xlab("Condition")

#pdf("Paper_1_Figures/Figure_3/jive_boxplots.pdf")
#print(pc1.box)
#print(pc2.box)
#dev.off()

p1 <- cowplot::plot_grid(pca.plot.points.connected.ellipse, pc2.box, pc1.box, 
                         align = "hv", rel_heights = c(1, .7), rel_widths = c(1, .4))

pdf(JPC.FIG.OUT.PATH, height = 10, width = 16)
print(p1)
dev.off()

mad_plot <- 
  ggplot(data = joint_sub, aes(color = cond_group)) +
  geom_abline(slope = 1, intercept = 0) + 
  geom_point(data = joint_sub, aes(x=mad_pc1, y=mad_pc2),size=5)+
  geom_text_repel(data = joint_sub, aes(x=mad_pc1, y=mad_pc2, label = cond.abbrev),size=5)+
  theme_bw() +
  xlab("MAD jPC1") +
  ylab("MAD jPC2") +
  xlim(c(0, 20)) + 
  ylim(c(0, 20)) + 
  scale_color_manual(values = scales::hue_pal()(4)[c(1,2,4)]) +
  geom_abline(slope = 1, intercept = 0) + 
  theme(legend.position = "none")
  #theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank()) +

ggsave(plot = mad_plot, filename = MAD.PLOT.OUT.PATH, height = 4.5, width = 4.5)


