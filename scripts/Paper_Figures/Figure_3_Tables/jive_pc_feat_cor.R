library(tidyverse)

if(!exists("snakemake")){
  setwd("../../..")
  source("scripts/util/paper/parse_snakemake.R")
  parse_snakemake(rule = "figure_3_jpc_feature_cor")
}

JIVE.PC.PATH <- snakemake@input[["jive_pca"]]
SOMA.PATH <- snakemake@input[["soma"]]
ARRAY.PATH <- snakemake@input[["array"]]


TAB.OUT.PATH <- snakemake@output[[1]]

prcomp.list <- readRDS(JIVE.PC.PATH)
soma.eset <- readRDS(SOMA.PATH)
array.eset <- readRDS(ARRAY.PATH)

#Get only the first three PC's of the joint. All other PC's essentially have eigen values of 0

joint <- prcomp.list$joint$x[, 1:3]
array_indiv = prcomp.list$array.ind$x[, 1:2]
soma_indiv = prcomp.list$soma.ind$x[, 1:2]

jive_pc_list <- list(joint = joint, transcriptome_individual = array_indiv, proteome_individual = soma_indiv)


#Put expressionsets into list and make sure that the patient_id are sampleNames/rownames of the expression matrix
eset.list <- list(protein = soma.eset, gene = array.eset)
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


get_cor_dat<- function(eset, jive_pc_mat, method){
  intersection <- intersect(rownames(jive_pc_mat), eset$patient_id)
  mat <- exprs(eset)
  mat <- mat[ ,match(intersection, eset$patient_id)]
  mat <- mat[complete.cases(mat),]
  mat <- t(mat)
  
  jive_pc_mat <- jive_pc_mat[match(intersection, rownames(jive_pc_mat)),]
  stopifnot(all.equal(rownames(mat), rownames(jive_pc_mat)))

  lapply(colnames(jive_pc_mat), function(PC){
    lapply(colnames(mat), function(feature){
      x <- jive_pc_mat[, PC]
      y <- mat[, feature]
      result <- do_cortest(x, y, method = method)
      data.frame(cor = result$estimate, p = result$p.value, PC = PC, feature = feature, stringsAsFactors = FALSE)
    }) %>% bind_rows()
  }) %>% bind_rows()
}

cordat_list <- lapply(jive_pc_list, function(jive_pc_mat){
   lapply(eset.list, get_cor_dat, jive_pc_mat = jive_pc_mat, method = "pearson")
})

cordat <- lapply(cordat_list, function(sublist){
  bind_rows(sublist, .id = "input_data")
}) %>% bind_rows(.id = "jive_pc_type")

cordat <- cordat %>%
        mutate(input_data = gsub("protein", "proteomics", input_data)) %>%
        mutate(input_data = gsub("gene", "transcriptomics", input_data))

cordat <- cordat %>%
        mutate(tmp = paste(jive_pc_type, input_data)) %>%
        filter(!tmp %in% c("transcriptome_individual proteomics", 
                           "proteome_individual transcriptomics")) #%>%
        #select(-tmp)

cordat <- cordat %>% select(-tmp)

cordat$PC[cordat$jive_pc_type == "joint"] <- 
        paste0("j", cordat$PC[cordat$jive_pc_type == "joint"])

cordat$PC[cordat$jive_pc_type == "transcriptome_individual"] <- 
        paste0("transcriptome_i", cordat$PC[cordat$jive_pc_type == "transcriptome_individual"])

cordat$PC[cordat$jive_pc_type == "proteome_individual"] <- 
        paste0("proteome_i", cordat$PC[cordat$jive_pc_type == "proteome_individual"])

write_csv(cordat, TAB.OUT.PATH)
