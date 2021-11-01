library(tidyverse)
library(gridExtra)


setwd("/hpcdata/sg/sg_data/PROJECTS/Monogenic_Project")

JIVE.PC.PATH <- "Integration_output/jive/subject/prcomp_list.rds"
TBNK.PATH <- "Data/TBNK/data_analysis_ready/tbnk_training_subject_level_eset.rds"
SOMA.PATH <- "Data/Somalogic/analysis_output/wgcna_results/scores_subject_level.rds"
ARRAY.PATH <- "Data/Microarray/analysis_output/WGCNA/array_subject_scores.rds"

prcomp.list <- readRDS(JIVE.PC.PATH)
tbnk.eset <- readRDS(TBNK.PATH)
soma.modules <- readRDS(SOMA.PATH)
array.modules <- readRDS(ARRAY.PATH)

#Get only the first three PC's of the joint. All other PC's essentially have eigen values of 0

array_indiv <- prcomp.list$array.ind$x[, 1:5]
soma_indiv <- prcomp.list$soma.ind$x[, 1:5]

#Put expressionsets into list and make sure that the patient_id are sampleNames/rownames of the expression matrix

eset.list <- list(soma = soma.modules, array = array.modules, tbnk = tbnk.eset)
eset.list <- lapply(eset.list, function(eset){
  sampleNames(eset) <- eset[["patient_id"]]
  eset
})

get_cormat <- function(eset, joint){
  intersection <- intersect(rownames(joint), eset$patient_id)
  mat <- exprs(eset)
  mat <- mat[ ,match(intersection, eset$patient_id)]
  mat <- mat[complete.cases(mat),]
  mat <- t(mat)
  
  joint <- joint[match(intersection, rownames(joint)),]
  stopifnot(all.equal(rownames(mat), rownames(joint)))
  
  cormat <- cor(mat, joint, method = "spearman")
  cormat
}

do_cortest <- function(x, y, method){
  intersection <- intersect(names(x), names(y))
  
  x <- x[match(intersection, names(x))]
  y <- y[match(intersection, names(y))]
  
  stopifnot(all.equal(names(x), names(y)))
  
  cor.test(x, y, method = method)
}

cormat_to_df <- function(cormat, main){
  cormat %>% 
  as.data.frame() %>% 
  mutate(feature = rownames(.)) %>%
  gather(key = PC, value = cor, -feature)
}

cormat.list.array <- lapply(eset.list, get_cormat, array_indiv)
cormat.list.soma <- lapply(eset.list, get_cormat, soma_indiv)

plot_dat_array <- lapply(cormat.list.array, cormat_to_df) %>%
        bind_rows(.id = "datatype")
plot_dat_soma <- lapply(cormat.list.soma, cormat_to_df) %>%
        bind_rows(.id = "datatype")

plot_dat <- list(array = plot_dat_array, soma = plot_dat_soma) %>%
        bind_rows(.id = "pc_source")


p <- ggplot(plot_dat, aes(x = PC, y = feature)) +
  geom_tile(aes(fill = cor)) +
  scale_radius(limits = c(0,1)) + 
  scale_fill_gradient2(low = "blue", mid = "white", 
                        high = "red", limits = c(-1, 1)) + 
  theme_classic() +
  facet_grid(datatype ~ pc_source, scales = "free", space = "free")
  #theme(axis.line = element_blank(), 
  #      axis.ticks = element_blank(), 
  #      axis.title = element_blank(),
  #      legend.position = "left") +

pdf("Paper_1_Figures/Supplemental_Figure_3/jive_individual_pc_cor.pdf", width = 8, height = 5)
print(p)
dev.off()


tmp <- readRDS("Paper_1_Figures/Figure_2/red_module_subcluster_enrich_top10.rds")
