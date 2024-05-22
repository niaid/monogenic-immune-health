library(tidyverse)
library(gridExtra)

#If there are issues with this script it is probably because I changed something with the featurename replacment


JIVE.PC.PATH <- snakemake@input[["jive_pcs"]]
TBNK.PATH <- snakemake@input[["tbnk"]]
SOMA.PATH <- snakemake@input[["soma_mod_scores"]]
ARRAY.PATH <- snakemake@input[["array_mod_scores"]]
#setwd("../../..")
#JIVE.PC.PATH <- "Pipeline_out/Integration_output/jive/subject/prcomp_list.rds"
#TBNK.PATH <- "Pipeline_out/Data/TBNK/data_analysis_ready/tbnk_training_subject_level_eset.rds"
#SOMA.PATH <- "Pipeline_out/Data/Somalogic/analysis_output/wgcna_results/scores_subject_level.rds"
#ARRAY.PATH <- "Pipeline_out/Data/Microarray/analysis_output/WGCNA/array_subject_scores.rds"


FIG.OUT.PATH <- snakemake@output[["figure"]]
TAB.OUT.PATH <- snakemake@output[["table"]]
#FIG.OUT.PATH <- "Pipeline_out/Paper_1_Figures/Figure_3/jive_pc_cor.pdf"

prcomp.list <- readRDS(JIVE.PC.PATH)
tbnk.eset <- readRDS(TBNK.PATH)
soma.modules <- readRDS(SOMA.PATH)
array.modules <- readRDS(ARRAY.PATH)

source("scripts/util/Plotting/tbnk_featurename_replace.R")
#Get only the first three PC's of the joint. All other PC's essentially have eigen values of 0

joint <- prcomp.list$joint$x[, 1:3]

#Put expressionsets into list and make sure that the patient_id are sampleNames/rownames of the expression matrix
eset.list <- list(protein = soma.modules, gene = array.modules, tbnk = tbnk.eset)
eset.list <- lapply(eset.list, function(eset){
  sampleNames(eset) <- eset[["patient_id"]]
  eset
})

do_cortest <- function(x, y, method){
  intersection <- intersect(names(x), names(y))
  
  x <- x[match(intersection, names(x))]
  y <- y[match(intersection, names(y))]
  
  stopifnot(all.equal(names(x), names(y)))
  
  cor.test(x, y, method = method)
}


get_cor_dat<- function(eset, joint, method){
  intersection <- intersect(rownames(joint), eset$patient_id)
  mat <- exprs(eset)
  mat <- mat[ ,match(intersection, eset$patient_id)]
  mat <- mat[complete.cases(mat),]
  mat <- t(mat)
  
  joint <- joint[match(intersection, rownames(joint)),]
  stopifnot(all.equal(rownames(mat), rownames(joint)))

  lapply(colnames(joint), function(PC){
    lapply(colnames(mat), function(feature){
      x <- joint[, PC]
      y <- mat[, feature]
      result <- do_cortest(x, y, method = method)
      data.frame(cor = result$estimate, p = result$p.value, PC = PC, feature = feature, stringsAsFactors = FALSE)
    }) %>% bind_rows()
  }) %>% bind_rows()
}

cordat_list <- lapply(eset.list, get_cor_dat, joint = joint, method = "spearman")

cordat_list[[1]] <- cordat_list[[1]] %>% mutate(feature2 = replace_mod_names_single_type(feature, "PM"))
cordat_list[[2]] <- cordat_list[[2]] %>% mutate(feature2 = replace_mod_names_single_type(feature, "TM"))
cordat_list[[3]] <- cordat_list[[3]] %>% mutate(feature2 = feature)

cordat <- bind_rows(cordat_list, .id = "feature_type")

cordat <- cordat %>%
        mutate(feature2 = replace(feature2, is.na(feature2), feature[is.na(feature2)])) %>% 
        mutate(p.adj = p.adjust(p, method = "fdr")) %>%
        mutate(asterisk = ifelse(p.adj < .05, "*", ""))

cordat$feature2 <- factor(cordat$feature2)

levels(cordat$feature2) <- replace_tbnk_names(levels(cordat$feature2))

cordat <- cordat %>% 
        mutate(feature_type2 = replace(feature_type, feature_type == "tbnk",
                                       tbnk_groups(feature2[feature_type == "tbnk"], 
                                                   "new name")))

lev_order <- c("TM", "PM", "Innate", "Lymphocytes", "RBC & PLT")
cordat <- cordat %>%
        mutate(feature_type2 = gsub("protein", "PM", feature_type2)) %>%
        mutate(feature_type2 = gsub("gene", "TM", feature_type2)) %>%
        mutate(feature_type2 = factor(feature_type2, levels = lev_order))

cordat <- cordat %>%
        mutate(PC = paste0("j", PC))


#write the table----
cordat %>% select(-c("asterisk", "feature_type", "feature")) %>%
        rename(feature_type = feature_type2, feature = feature2) %>%
        write_csv(TAB.OUT.PATH)

#plot the figure ---
p <- ggplot(cordat %>% filter(PC %in% c("jPC1", "jPC2")), aes(y = PC, x = feature2)) + geom_tile(aes(fill = cor)) +
  scale_radius(limits = c(0,1)) + 
  scale_fill_gradient2(low = "blue", mid = "white", 
                        high = "red", limits = c(-1, 1)) + 
  geom_text(aes(label = asterisk), color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.background.y = element_blank(), strip.text.y = element_blank(),
        axis.line = element_blank(), 
        #strip.background.x = element_rect(colour="black", fill="grey90"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_grid(1 ~ feature_type2, scales = "free", space = "free") +
  ylab("") + xlab("")
  #theme(axis.line = element_blank(), 
  #      axis.ticks = element_blank(), 
  #      axis.title = element_blank(),
  #      legend.position = "left") +

pdf(FIG.OUT.PATH, width = 7, height = 3.7)
print(p)
dev.off()


