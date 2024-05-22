library(tidyverse)
library(Biobase)
library(ggpubr)

if(!exists("snakemake")){
  setwd("../../..")
  source("scripts/util/paper/parse_snakemake.R")
  parse_snakemake(rule = "figure_2_il23_boxplot")
}
#-------
#setwd("../../..")
#PROTEIN_IN_PATH <- "Pipeline_out/Data/Somalogic/data_analysis_ready/analysis_ready_subject_level_training_somalogic.rds"
#SOMALOGIC_DE_RESULTS_IN_PATH <- 'Pipeline_out/Data/Somalogic/analysis_output/differential_expression/somalogic_features_DE_results.rds'
#FIG.OUT.PATH <- "Pipeline_out/Paper_1_Figures/Figure_2/IL23.pdf"
#-------

PROTEIN_IN_PATH <- snakemake@input[["soma_feat_data"]]#"Data/Somalogic/data_analysis_ready/analysis_ready_subject_level_training_somalogic.rds"

SOMALOGIC_DE_RESULTS_IN_PATH <- snakemake@input[["soma_feat_de"]]#'Data/Somalogic/analysis_output/differential_expression/somalogic_features_DE_results.rds'

FIG.OUT.PATH <- snakemake@output[[1]]#"Paper_1_Figures/Figure_2/IL23.pdf"

source("scripts/util/paper/abbrev_cond.R")

protein_results <- readRDS(SOMALOGIC_DE_RESULTS_IN_PATH)

protein_eset <- readRDS(PROTEIN_IN_PATH)

soma_q_dat <- protein_results$versus.healthy$adj.P.Val %>%
        as.data.frame() %>%
        rownames_to_column(var = "feature") %>%
        gather(key = "condition", value = "q_val", -feature)

asterisk_vec <- function(p_vec){
  out_vec <- rep("", length(p_vec))
  out_vec <- replace(out_vec, p_vec < .05, "*")
  out_vec <- replace(out_vec, p_vec < .01, "**")
  out_vec <- replace(out_vec, p_vec < .001, "***")
  out_vec
}

IL23_signif <- soma_q_dat %>%
        filter(feature =="IL.23") %>%
        mutate(IL23_signif= asterisk_vec(q_val))


protein_meta <- protein_eset %>%
        pData() %>%
        select(patient_id, gender, Age, condition) %>%
        mutate(IFN.g.protein = exprs(protein_eset)["IFN.g", ]) %>%
        mutate(IL23 = exprs(protein_eset)["IL.23", ])

protein_meta <- protein_meta %>%
        left_join(IL23_signif %>% select(condition, IL23_signif))

plot_dat <- protein_meta

IL23_medians <- plot_dat %>% group_by(condition) %>%
        summarise(med_IL23 = median(IL23)) %>%
        deframe()

plot_dat <- plot_dat %>%
        mutate(cond_group = group_cond(condition))

plot_dat <- plot_dat %>%
        mutate(condition = factor(condition, levels = names(sort(IL23_medians))))


levels(plot_dat$condition) <- abbrev_cond(levels(plot_dat$condition))

pdf(FIG.OUT.PATH, height = 3.5, width = 3.5)
#p <- ggplot(plot_dat, aes(x = condition, y = IL23, fill = cond_group)) +
#        geom_boxplot(outlier.shape = NA) +
#        #geom_jitter(height = 0) +
#        ggbeeswarm::geom_beeswarm() +
#        geom_text(aes(x = condition, y = max(IL23) *1.1, label = IL23_signif), size = 6, color = "red", position = position_nudge(x = -.5)) +
#        ylim(min(plot_dat$IL23), max(plot_dat$IL23) * 1.12) +
#        facet_grid(condition == "Healthy" ~ 1, scales = "free_y", space = "free_y") +
#        coord_flip() +
#        theme_bw() +
#        theme(legend.position = "none", 
#              strip.background = element_blank(),
#              strip.text.x = element_blank(), 
#              strip.text.y = element_blank(),
#              panel.spacing = unit(0, "lines"))
#print(p)


p <- ggplot(plot_dat, aes(y = IL23, x = reorder(condition, IL23, median))) +
                               #color = condition == "Healthy")) +
        #geom_boxplot(outlier.shape = NA) +
        geom_boxplot(outlier.shape = NA, aes(fill = cond_group)) +
        geom_text(angle=90, aes(x = condition, y = max(IL23) * 1.1, label = IL23_signif), size = 6, color = "red", position = position_nudge(x = -.5)) +
        ylim(min(plot_dat$IL23),max(plot_dat$IL23) * 1.2) +
        #geom_jitter(width = 0, alpha = .5) +
        ggbeeswarm::geom_beeswarm(alpha = .5, size = 1) +
        xlab("IL23") +
        facet_grid(1 ~condition == "Healthy", scales = "free_x", space = "free_x") +
        theme_bw() +
        ylab("") + 
        xlab("") + 
        #scale_color_manual(values = c("black", "red")) + 
        theme(legend.position = "none", 
              strip.background = element_blank(),
              strip.text.x = element_blank(), 
              strip.text.y = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              panel.spacing = unit(0, "lines"))

print(p)
dev.off()


