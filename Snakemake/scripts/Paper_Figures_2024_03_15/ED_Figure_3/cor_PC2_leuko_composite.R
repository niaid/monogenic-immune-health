suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
  library(tidyverse)
  library(ggpubr)
})



source("scripts/util/Processing/averageRepeatSamples.R")
source("scripts/util/paper/abbrev_cond.R")

if(!exists("snakemake")){
  setwd("../../..")
  source("scripts/util/paper/parse_snakemake.R")
  parse_snakemake(rule = "supplemental_figure_3_pc2_leuko_composite")
}

JIVE.PATH <- snakemake@input[["jive"]]#"Integration_output/jive/subject/jive.rds"
JIVE.PC.PATH <- snakemake@input[["jive_pcs"]]#"Integration_output/jive/subject/prcomp_list.rds"
TBNK.PATH <- snakemake@input[["tbnk"]]#"Data/TBNK/data_analysis_ready/tbnk_training_subject_level_eset.rds"

PLOT_OUT_PATH <- snakemake@output[[1]]#"Paper_1_Figures/Supplemental_Figure_3/leuko_composite_pc2_cor.pdf"
MARROW_PLOT_OUT_PATH <- snakemake@output[[2]]#"Paper_1_Figures/Supplemental_Figure_3/leuko_composite_pc2_cor.pdf"
PLOT_SEPARATE_OUT_PATH <- snakemake@output[[3]]#"Paper_1_Figures/Supplemental_Figure_3/leuko_composite_separate_pc2_cor.pdf"

tbnk.eset <- readRDS(TBNK.PATH)

prcomp.list <- readRDS(JIVE.PC.PATH)
joint <- prcomp.list$joint$x

jive <- readRDS(JIVE.PATH)
pdat <- jive$pdat

intersection <- intersect(rownames(joint), tbnk.eset$patient_id)
tbnk.mat <- exprs(tbnk.eset)
tbnk.mat <- tbnk.mat[,match(intersection, tbnk.eset$patient_id)]
tbnk.mat <- t(tbnk.mat)

joint <- joint[match(intersection, rownames(joint)),]

stopifnot(all.equal(rownames(tbnk.mat), rownames(joint)))
pdat <- pdat %>% filter(patient_id %in% rownames(joint))

#Select lymphocytes, monocytes, neutrophils absolute counts

keep.cells <- c("neutrophil_abs", "monocytes_abs", "lymphocytes_abs")
tbnk.mat <- tbnk.mat[, keep.cells]

# Make everything z score for healthy mean and sd

healthy.means <- apply(tbnk.mat[pdat$condition == "Healthy",], 2, mean)
healthy.sd <- apply(tbnk.mat[pdat$condition == "Healthy",], 2, sd)

tbnk.z <- tbnk.mat
for(i in seq_len(ncol(tbnk.z))){
  tbnk.z[, i] <- (tbnk.z[, i] - healthy.means[[i]]) / healthy.means[[i]]
}

#Create composite score

#Average of the Z-scores
tbnk.composite <- apply(tbnk.z, 1, mean)

dat <- pdat %>%
  mutate(composite.score = tbnk.composite) %>%
  mutate(PC2 = joint[, "PC2"]) %>%
  mutate(cond.grouped = group_cond(condition)) %>%
  mutate(cond.abbrev = abbrev_cond(condition))

#dat %>% 
#  group_by(condition) %>%
#  summarise(pc2.med = median(PC2)) %>% 
#  arrange(pc2.med)

#These are the conditions that will be included in the scatter of PC2 vs composite score when show each condition in facets


conditions.of.interest <- c("Healthy", "DADA2", "GATA2", "CTLA4", "PGM3", "PI3K", "TERC", "TERT")


#Add column that can be used to select conditions of interest and add annotation that groups the Terts and Tercs


dat <- 
  dat %>%
  mutate(condition2 = replace(condition, which(!condition %in% conditions.of.interest), "other")) %>% 
  mutate(condition2 = replace(condition2, condition %in% c("TERT", "TERC"), "TERT/TERC"))


#Plot across everyone

pdf(PLOT_OUT_PATH, height =5, width = 5)
ggplot(dat, aes(x = composite.score, y = PC2)) + 
  geom_text(aes(color = cond.abbrev, label = cond.abbrev), size = 2) + 
  ylab("jPC2") +
  stat_cor(method = "spearman") +
  theme_bw() +
  theme(legend.position = "none")
dev.off()


#2nd plot. just gata2
normal_pats <- c("P129", "P164", "P182", "P150")
mild_pats <- c("P97", "P101")
mds_pats <- c("P168", "P86", "P166")

gata2_dat <- dat %>% 
        filter(condition == "GATA2") %>%
        mutate(marrow_status = NA) %>%
        mutate(marrow_status = replace(marrow_status, patient_id %in% normal_pats, "normal")) %>%
        mutate(marrow_status = replace(marrow_status, patient_id %in% mild_pats, "mild G2BMD")) %>%
        mutate(marrow_status = replace(marrow_status, patient_id %in% mds_pats, "MDS")) %>%
        mutate(marrow_status = factor(marrow_status, levels = c("normal", "mild G2BMD", "MDS")))

p <- ggplot(gata2_dat, aes(x = composite.score, y = PC2)) +
        geom_point(aes(shape = marrow_status, color = marrow_status)) +
        scale_color_manual(values = c("black", "orange", "red")) +
        labs(color = "Marrow Status", shape = "Marrow Status") +
        theme_bw() +
        facet_wrap(~"GATA2")

ggsave(plot = p, filename = MARROW_PLOT_OUT_PATH, height = 2, width = 3.5)

#Plot by condition with p values
p <- dat %>%
  filter(condition %in% conditions.of.interest) %>%
  ggplot(aes(x = composite.score, y = PC2)) + 
  geom_point(aes(color = cond.abbrev)) + 
  ylab("jPC2") +
  stat_cor(method = "spearman") + 
  facet_wrap(~condition2, nrow = 4) + 
  theme_bw() +
  theme(legend.position = "none")
ggsave(plot = p, filename = PLOT_SEPARATE_OUT_PATH, width = 4)


