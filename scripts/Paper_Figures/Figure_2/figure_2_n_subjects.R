library(tidyverse)
library(Biobase)
library(BiocGenerics)

ESETS.IN.PATHS = list(
  somalogic.modules = snakemake@input[[1]],# 'Data/Somalogic/analysis_output/wgcna_results/scores_subject_level.rds',
  somalogic.features = snakemake@input[[2]],# 'Data/Somalogic/data_analysis_ready/analysis_ready_subject_level_training_somalogic.rds',
  microarray.modules = snakemake@input[[3]],# 'Data/Microarray/analysis_output/WGCNA/array_subject_scores.rds',
  microarray.features = snakemake@input[[4]],# 'Data/Microarray/data_analysis_ready/eset_batch_training_subject.rds',
  tbnks = snakemake@input[[5]]# 'Data/TBNK/data_analysis_ready/tbnk_training_subject_level_eset.rds'
)
#ESETS.IN.PATHS = list(
#  somalogic.modules =  'Data/Somalogic/analysis_output/wgcna_results/scores_subject_level.rds',
#  somalogic.features = 'Data/Somalogic/data_analysis_ready/analysis_ready_subject_level_training_somalogic.rds',
#  microarray.modules = 'Data/Microarray/analysis_output/WGCNA/array_subject_scores.rds',
#  microarray.features = 'Data/Microarray/data_analysis_ready/eset_batch_training_subject.rds',
#  tbnks = 'Data/TBNK/data_analysis_ready/tbnk_training_subject_level_eset.rds'
#)

source("scripts/util/paper/abbrev_cond.R")

FIG.OUT.PATH <- snakemake@output[[1]]#"Paper_1_Figures/Figure_2/all_feature_bubble_plot_condition_separated_nsubj_per_comparison.pdf"

esets <- lapply(ESETS.IN.PATHS, readRDS)

# Get the patient-condition pairings from each eset in order to calculate the number of subjects coming from each data type
sample.datas = lapply(esets, function(eset) {
  sample.data = pData(eset)
  sample.data[, c('patient_id','condition')]
})

n_subj_dat <- lapply(sample.datas, function(dat){
  dat %>% group_by(condition) %>%
          summarise(n = n())
}) %>% bind_rows(.id = "feat_category")

n_subj_dat <- n_subj_dat %>%
        filter(feat_category %in% c("somalogic.modules", "microarray.modules", "tbnks")) %>%
        mutate(feat_category = gsub("somalogic.modules", "protein", feat_category)) %>%
        mutate(feat_category = gsub("microarray.modules", "gene", feat_category)) %>%
        mutate(feat_category = gsub("protein", "PM", feat_category)) %>%
        mutate(feat_category = gsub("gene", "TM", feat_category)) %>%
        mutate(feat_category = gsub("tbnks", "CBC + TBNK", feat_category))

sum_n_subj <- n_subj_dat %>%
        group_by(condition) %>%
        summarise(sum_n_subj = sum(n)) %>%
        arrange(sum_n_subj)

cond_order <- sum_n_subj$condition

n_subj_dat <- n_subj_dat %>% mutate(cond_group = group_cond(condition)) %>%
        left_join(sum_n_subj) %>%
        mutate(condition = factor(condition, levels = cond_order))

levels(n_subj_dat$condition) <- abbrev_cond(levels(n_subj_dat$condition))

pdf(FIG.OUT.PATH, height = 2.5, width = 6)
p2 <- ggplot(n_subj_dat %>% filter( !condition %in% c("TERT.TERC", "AI", "PID")), 
             aes(x = condition, feat_category)) +
        geom_tile(aes(fill= log2(n))) +
        geom_text(aes(label = n), color = "white", angle = 90) + 
        facet_grid(1 ~ cond_group, scales = "free", space = "free") +
        ggtitle("number of subjects in each comparison") +
        #scale_fill_viridis_c()+
        scale_fill_gradient(low = "white", high = "black") +
        theme_bw() +
        ylab("") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              strip.background = element_blank(), strip.text = element_blank())
print(p2)
dev.off()
