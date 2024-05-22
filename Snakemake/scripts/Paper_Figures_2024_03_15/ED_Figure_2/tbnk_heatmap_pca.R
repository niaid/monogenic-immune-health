library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggfortify)
library(ggrepel)
library(Biobase)

source("scripts/util/paper/abbrev_cond.R")
source("scripts/util/Plotting/tbnk_featurename_replace.R")

TBNK.ESET.IN.PATH <- snakemake@input[[1]]#"Data/TBNK/data_analysis_ready/tbnk_training_subject_level_eset.rds"
eset <- readRDS(TBNK.ESET.IN.PATH)

HEATMAP.OUT.PATH <- snakemake@output[["heatmap"]]#"Paper_1_Figures/Supplemental_Figure_2/tbnk_heatmap.pdf"
PCA.OUT.PATH <- snakemake@output[["pca"]]#"Paper_1_Figures/Supplemental_Figure_2/tbnk_pca.pdf"
#keep.features <- c(
#  "cd4_cd3_count","nk_cells_count","cd3_count","cd8_cd3_count",
#  "cd19_count","wbc","rbc",
#  "hemoglobin","mcv",
#  "mch","rdw",
#  "platelet_count","neutrophil_abs",
#  "lymphocytes_abs","monocytes_abs","eosinophil_abs",
#  "basophil_abs"
#)

mat <- exprs(eset[1:18, ]) %>% t()
pr_obj <- prcomp(mat, scale = TRUE, center = TRUE)


meta <- pData(eset) %>%
  mutate(cond.abbrev = abbrev_cond(condition)) %>%
  mutate(cond.grouped = group_cond(condition)) %>%
  mutate(PC1 = pr_obj$x[, "PC1"],
         PC2 = pr_obj$x[, "PC2"])

centroids <- meta %>% group_by(condition) %>%
        summarise(mean_pc1 = mean(PC1), mean_pc2 = mean(PC2),
        n_subj = n())

meta <- left_join(meta, centroids)

centroid_sub <- meta %>% 
        filter(n_subj > 3) %>%
        select(mean_pc1, mean_pc2, cond.grouped, cond.abbrev) %>% distinct()

pca.plot.points <- 
  ggplot(meta, aes(x = PC1, y = PC2, color = cond.grouped)) +
  geom_text(aes(label = cond.abbrev), size = 2) +
  geom_point(data = centroid_sub, aes(x=mean_pc1, y=mean_pc2),size=5)+
  geom_segment(aes(x=mean_pc1, y=mean_pc2, xend=PC1, yend=PC2), alpha = .2) +
  geom_text_repel(data = centroid_sub, aes(x=mean_pc1, y=mean_pc2, label = cond.abbrev),size=5)+
  theme_bw() +
  theme(legend.position = "none")

#pca.plot.age <- 
#  #ggplot(meta, aes(x = PC1, y = PC2, color = cond.abbrev)) +
#  ggplot(meta, aes(x = PC1, y = PC2)) +
#  geom_text(aes(label = cond.abbrev, color = Age), size = 2) +
#  #geom_point(aes(x=mean_pc1, y=mean_pc2),size=5)+
#  #geom_segment(aes(x=mean_pc1, y=mean_pc2, xend=PC1, yend=PC2), alpha = .2) +
#  #geom_text_repel(data = centroid_sub, aes(x=mean_pc1, y=mean_pc2, label = cond.abbrev),size=5)+
#  theme_bw()

#file.remove("Paper_1_Figures/Supplemental_Figure_2/tbnk_pca.pdf")
pdf(PCA.OUT.PATH, height = 5, width = 5.5)
print(pca.plot.points)
#print(pca.plot.age)
dev.off()

#Heatmap showing Z-score of CBC parameters

#The big blue group clustering together is Job; they have high eosinophils


condition <- abbrev_cond(eset$condition)
conditions <- table(condition)
large.conditions <- conditions[table(condition) > 10]
large.condition <- condition
large.condition[!large.condition %in% names(large.conditions)] <- "Other"
large.condition <- factor(large.condition)

#All subjects

#cond <- factor(eset$condition)
annotation <- data.frame(All_groups = condition, Major_groups = large.condition, age = eset$Age)
rownames(annotation) <- colnames(eset)

breaksList = seq(-3, 3, by = .01)
pdf(HEATMAP.OUT.PATH, height = 9, width = 14)
pheatmap(exprs(eset) %>% t %>% scale %>% t %>% `rownames<-`(replace_tbnk_names(rownames(.))),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList,
         show_colnames = FALSE, 
         annotation_col = annotation)
dev.off()

