library(aricode)
library(tidyverse)

source("jaccard_matrix.R")

soma_modules_official <- readRDS("../../Snakemake/Pipeline_out/Data/Somalogic/analysis_output/wgcna_results/modules.rds")
array_modules_official <- readRDS("../../Snakemake/Pipeline_out/Data/Microarray/analysis_output/WGCNA/modules.rds")

soma_path <- "snakemake_workflow/varpart_objects/soma/"
soma_files <- list.files(soma_path, pattern = "modules")

array_path <- "snakemake_workflow/varpart_objects/array/"
array_files <- list.files(array_path, pattern = "modules")

soma_modules <- lapply(soma_files, function(x){readRDS(file.path(soma_path, x))})
array_modules <- lapply(array_files, function(x){readRDS(file.path(array_path, x))})


soma_ari <- sapply(soma_modules, function(x){ARI(x, soma_modules_official)})
soma_ami <- sapply(soma_modules, function(x){AMI(x, soma_modules_official)})

array_ari <- sapply(array_modules, function(x){ARI(x, array_modules_official)})
array_ami <- sapply(array_modules, function(x){AMI(x, array_modules_official)})

dat_list <- list(
  data.frame(data.type = "Protein", measure = "Adjusted Rand Index", value = soma_ari),
  data.frame(data.type = "Protein", measure = "Adjusted Mutual Information", value = soma_ami),
  data.frame(data.type = "Gene Expression", measure = "Adjusted Rand Index", value = array_ari),
  data.frame(data.type = "Gene Expression", measure = "Adjusted Mutual Information", value = array_ami)
)

dat <- bind_rows(dat_list)

p <- ggplot(dat, aes(x = data.type, y = value)) +
        facet_wrap(~measure, nrow=1) +
        geom_boxplot() +
        ylim(0,1) +
        theme_bw()
dir.create("plots")
ggsave(plot = p, filename = "plots/rand_index_and_mutual_info_modules_jackknife.pdf", height= 4)

convert_to_list <- function(mod_vector){
  all_mods <- unique(mod_vector)
  names(all_mods) <- all_mods

  lapply(all_mods, function(mod){
    names(mod_vector)[mod_vector == mod]
  })
}

soma_modules_official_list <- convert_to_list(soma_modules_official)
array_modules_official_list <- convert_to_list(array_modules_official)

soma_modules_list <- lapply(soma_modules, convert_to_list)
array_modules_list <- lapply(array_modules, convert_to_list)

soma_jaccard_list <- lapply(soma_modules_list, function(mod_list){
  jacc_mat <- jaccardMat2(soma_modules_official_list, mod_list)
  apply(jacc_mat, 1, max)
})
soma_jaccard_dat <- lapply(soma_jaccard_list, function(x){
  enframe(x, name = "module_color", value = "jaccard_max")
}) %>% bind_rows()

array_jaccard_list <- lapply(array_modules_list, function(mod_list){
  jacc_mat <- jaccardMat2(array_modules_official_list, mod_list)
  apply(jacc_mat, 1, max)
})
array_jaccard_dat <- lapply(array_jaccard_list, function(x){
  enframe(x, name = "module_color", value = "jaccard_max")
}) %>% bind_rows()

source("../../Snakemake/scripts/util/Plotting/tbnk_featurename_replace.R")

array_jaccard_dat <- array_jaccard_dat %>%
        mutate(module = 
               replace_mod_names_single_type(module_color, 1, excelpath = "../../Snakemake/Inputs/Reference/module_name_map.xlsx")
       )

soma_jaccard_dat <- soma_jaccard_dat %>%
        mutate(module = 
               replace_mod_names_single_type(module_color, 2, excelpath = "../../Snakemake/Inputs/Reference/module_name_map.xlsx")
       )

pdf("plots/max_jaccard_overlap_modules_jackknife.pdf", height = 4, width =4)
p <- ggplot(array_jaccard_dat, aes(x = module, y = jaccard_max)) +
        geom_boxplot() +
        ylim(c(0,1)) +
        theme_bw() +
        ggtitle("TMs")+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p)

p <- ggplot(soma_jaccard_dat, aes(x = module, y = jaccard_max)) +
        geom_boxplot() +
        ylim(c(0,1)) +
        theme_bw() +
        ggtitle("PMs") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p)
dev.off()

