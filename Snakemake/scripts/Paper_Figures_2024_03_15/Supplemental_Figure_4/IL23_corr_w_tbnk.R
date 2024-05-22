library(tidyverse)
library(Biobase)
library(ggpubr)

if(!exists("snakemake")){
  setwd("../../..")
  source("scripts/util/paper/parse_snakemake.R")
  parse_snakemake(rule = "figure_2_il23_tbnk_cor_dada2")
}


TBNK_IN_PATH <- snakemake@input[["tbnk_data"]]#"Data/TBNK/data_analysis_ready/tbnk_training_sample_level_eset.rds"

PROTEIN_IN_PATH <- snakemake@input[["soma_feat_data"]]#"Data/Somalogic/data_analysis_ready/analysis_ready_subject_level_training_somalogic.rds"

FIG.OUT.PATH <- snakemake@output[[1]]#"Paper_1_Figures/Figure_2/IL23_tbnk_cor_select_features.pdf"

tbnk_eset <- readRDS(TBNK_IN_PATH)
protein_eset <- readRDS(PROTEIN_IN_PATH)

tbnk_meta <- tbnk_eset %>%
        pData() %>%
        select(patient_id, gender, condition)

keep_tbnk_features <- grep("abs", featureNames(tbnk_eset), value = TRUE)
tbnk_meta <- bind_cols(tbnk_meta, as.data.frame(t(exprs(tbnk_eset)[keep_tbnk_features, ])))

protein_meta <- protein_eset %>%
        pData() %>%
        #select(patient_id, gender, Age, condition) %>%
        select(patient_id, gender, condition) %>%
        mutate(IFN.g.protein = exprs(protein_eset)["IFN.g", ]) %>%
        mutate(IL23 = exprs(protein_eset)["IL.23", ])

protein_meta %>%
        filter(condition == "DADA2") %>% 
        arrange(IL23)


plot_dat <- left_join(protein_meta, tbnk_meta)

plot_dat_long <- plot_dat %>%
        gather(key = "feature", value = "value", -c("patient_id", "gender", "condition", "IL23"))
        #gather(key = "feature", value = "value", -c("patient_id", "gender", "Age", "condition", "IL23"))

keep_features <- c("platelet_abs", "neutrophil_abs", "cd19_abs", "IFN.g.protein")

#df$genes <- factor(df$genes, levels = c("BA","MLL","pos","neg","PMLalpha+"),
# ordered = TRUE, labels=c("BA","MLL","pos","neg",expression(paste("PML", alpha,"+"))))
labs <-  c(expression(paste("PLT")), 
           expression(paste("Neutrophil (#)")), 
           expression(paste("CD19+ B Cells (#)")), 
           expression(paste("IFN-", gamma, " Protein")))
#labs <- c(lev[1:3], )
plot_dat_long_sub <- plot_dat_long %>% 
        filter(condition == "DADA2", feature %in% keep_features)%>%
        mutate(feature = factor(feature, levels = keep_features, labels = labs, ordered = TRUE))

pdf(FIG.OUT.PATH, height = 4, width = 4)
p2 <- ggplot(plot_dat_long_sub, aes(x = IL23, y = value)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE) +
        stat_cor(size = 3) +
        facet_wrap(~ feature, scales = "free", nrow = 2, labeller = label_parsed) +
        theme_bw() +
        ggtitle("DADA2 only")
print(p2)
dev.off()
