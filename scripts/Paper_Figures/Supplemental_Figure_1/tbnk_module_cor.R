library(tidyverse)
library(gridExtra)


#setwd("/hpcdata/sg/sg_data/PROJECTS/Monogenic_Project")

TBNK.PATH <- snakemake@input[["tbnk"]]#"Data/TBNK/data_analysis_ready/tbnk_training_subject_level_eset.rds"
SOMA.PATH <- snakemake@input[["soma"]]#"Data/Somalogic/analysis_output/wgcna_results/scores_subject_level.rds"
ARRAY.PATH <- snakemake@input[["array"]]#"Data/Microarray/analysis_output/WGCNA/array_subject_scores.rds"

FIG.OUT.PATH <- snakemake@output[[1]]#"Paper_1_Figures/Supplemental_Figure_1/tbnk_cor_w_modules2.pdf"

tbnk.eset <- readRDS(TBNK.PATH)
soma.modules <- readRDS(SOMA.PATH)
array.modules <- readRDS(ARRAY.PATH)

#rename modules
source("scripts/util/Plotting/tbnk_featurename_replace.R")

#Put expressionsets into list and make sure that the patient_id are sampleNames/rownames of the expression matrix

eset.list <- list(protein = soma.modules, gene = array.modules, tbnk = tbnk.eset)
eset.list <- lapply(eset.list, function(eset){
  sampleNames(eset) <- eset[["patient_id"]]
  eset
})
featureNames(eset.list$protein) <- replace_mod_names_single_type(featureNames(eset.list$protein), "PM")

featureNames(eset.list$gene) <- replace_mod_names_single_type(featureNames(eset.list$gene), "TM")


tbnk_mat <- exprs(eset.list$tbnk) %>% t()

do_cortest <- function(x, y, method){
  intersection <- intersect(names(x), names(y))
  
  x <- x[match(intersection, names(x))]
  y <- y[match(intersection, names(y))]
  
  stopifnot(all.equal(names(x), names(y)))
  
  cor.test(x, y, method = method)
}


get_cor_dat<- function(eset, mat2, method){
  intersection <- intersect(rownames(mat2), eset$patient_id)
  mat <- exprs(eset)
  mat <- mat[ ,match(intersection, eset$patient_id)]
  mat <- mat[complete.cases(mat),]
  mat <- t(mat)
  
  mat2 <- mat2[match(intersection, rownames(mat2)),]
  stopifnot(all.equal(rownames(mat), rownames(mat2)))

  lapply(colnames(mat2), function(feature_x){
    lapply(colnames(mat), function(feature_y){
      x <- mat2[, feature_x]
      y <- mat[, feature_y]
      result <- do_cortest(x, y, method = method)
      data.frame(cor = result$estimate, p = result$p.value, feature_x = feature_x, feature_y = feature_y)
    }) %>% bind_rows()
  }) %>% bind_rows()
}

cordat.list <- lapply(eset.list[c(1,2)], get_cor_dat, tbnk_mat, method = "spearman")

reorder_levels <- function(dat){
  cormat <- dat %>% 
          select(cor, feature_x, feature_y) %>%
          spread(key = feature_y, value = cor) %>%
          `rownames<-`(.$feature_x) %>%
          select(-feature_x) %>%
          as.matrix()

  hc_row <- cormat %>%
          dist() %>% hclust()
  dat$feature_x <- factor(dat$feature_x, levels = hc_row$labels[hc_row$order])

  hc_col <- cormat %>%
          t() %>%
          dist() %>% hclust()
  dat$feature_y <- factor(dat$feature_y, levels = hc_col$labels[hc_col$order])

  dat
}

cordat.list <- lapply(cordat.list, reorder_levels)

cordat.list <- lapply(cordat.list, function(dat){
  levels(dat$feature_x) <- replace_tbnk_names(levels(dat$feature_x))
  dat$feat_group <- tbnk_groups(dat$feature_x, "new name")

  dat
})

add_signif <- function(dat, method = "fdr", cutoff = .05){
  dat <- dat %>%
          mutate(p.adj = p.adjust(p, method = "fdr")) %>%
          mutate(asterisk = ifelse(p.adj < cutoff, "*", ""))
  dat
}

cordat.list <- lapply(cordat.list, add_signif)

pdf(FIG.OUT.PATH, height = 6, width = 4)
p <- ggplot(cordat.list[[1]] , aes(x = feature_y, y = feature_x)) + geom_tile(aes(fill = cor)) +
  #scale_radius(limits = c(0,1)) + 
  scale_fill_gradient2(low = "blue", mid = "white", 
                        high = "red", limits = c(-1, 1)) + 
  geom_text(aes(label = asterisk), color = "black") +
  facet_grid(feat_group~"PM", scales = "free_y", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.line = element_blank()) +
        #strip.text.x = element_blank(), 
        #strip.background.x = element_blank()) +
  xlab("") + ylab("")
print(p)

p <- ggplot(cordat.list[[2]] , aes(x = feature_y, y = feature_x)) + geom_tile(aes(fill = cor)) +
  #scale_radius(limits = c(0,1)) + 
  scale_fill_gradient2(low = "blue", mid = "white", 
                        high = "red", limits = c(-1, 1)) + 
  geom_text(aes(label = asterisk), color = "black") +
  theme_bw() +
  facet_grid(feat_group~"TM", scales = "free_y", space = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.line = element_blank()) + 
        #strip.text.x = element_blank(), 
        #strip.background.x = element_blank()) +
  xlab("") + ylab("")

print(p)
dev.off()

