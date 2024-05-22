library(tidyverse)

soma_modules_official <- readRDS("../../Snakemake/Pipeline_out/Data/Somalogic/analysis_output/wgcna_results/scores_subject_level.rds")
array_modules_official <- readRDS("../../Snakemake/Snakemake/Pipeline_out/Data/Microarray/analysis_output/WGCNA/array_subject_scores.rds")

soma_path <- "snakemake_workflow/varpart_objects/soma/"
soma_files <- list.files(soma_path, pattern = "subject_level.rds")

array_path <- "snakemake_workflow/varpart_objects/array/"
array_files <- list.files(array_path, pattern = "subject_level.rds")

soma_modules <- lapply(soma_files, function(x){readRDS(file.path(soma_path, x))})
array_modules <- lapply(array_files, function(x){readRDS(file.path(array_path, x))})

expr_soma_official <- t(exprs(soma_modules_official))
expr_array_official <- t(exprs(array_modules_official))

soma_cor_list <- lapply(soma_modules, function(modules_eset){
  mat <- t(exprs(modules_eset))
  cor_mat <- cor(expr_soma_official, mat)
  apply(cor_mat, 1, max)
})
soma_cor_dat <- lapply(soma_cor_list, function(x){
  enframe(x, name = "module_color", value = "cor_max")
}) %>% bind_rows()

array_cor_list <- lapply(array_modules, function(modules_eset){
  mat <- t(exprs(modules_eset))
  cor_mat <- cor(expr_array_official, mat)
  apply(cor_mat, 1, max)
})
array_cor_dat <- lapply(array_cor_list, function(x){
  enframe(x, name = "module_color", value = "cor_max")
}) %>% bind_rows()

source("../../Snakemake/scripts/util/Plotting/tbnk_featurename_replace.R")

array_cor_dat <- array_cor_dat %>%
        mutate(module = 
               replace_mod_names_single_type(module_color, 1, excelpath = "../../Snakemake/Inputs/Reference/module_name_map.xlsx")
       )

soma_cor_dat <- soma_cor_dat %>%
        mutate(module = 
               replace_mod_names_single_type(module_color, 2, excelpath = "../../Snakemake/Inputs/Reference/module_name_map.xlsx")
       )

pdf("plots/max_cor_overlap_modules_jackknife.pdf", height = 4, width =4)
p <- ggplot(array_cor_dat, aes(x = module, y = cor_max)) +
        geom_boxplot() +
        ylim(c(0,1)) +
        theme_bw() +
        ggtitle("TMs")+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p)

p <- ggplot(soma_cor_dat, aes(x = module, y = cor_max)) +
        geom_boxplot() +
        ylim(c(0,1)) +
        theme_bw() +
        ggtitle("PMs") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p)
dev.off()

