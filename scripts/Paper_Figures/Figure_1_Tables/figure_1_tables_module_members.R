library(tidyverse)
library(Biobase)
library(BiocGenerics)

#setwd("/hpcdata/sg/sg_data/PROJECTS/Monogenic_Project/")

if(!exists("snakemake")){
  setwd("../../..")
  source("scripts/util/paper/parse_snakemake.R")
  parse_snakemake("figure_1_tables_module_members")
}

SOMALOGIC_MODULE_IN_PATH <- snakemake@input[["soma_mod"]]#"Data/Somalogic/analysis_output/wgcna_results/modules.rds"
MICROARRAY_MODULE_IN_PATH <- snakemake@input[["array_mod"]]#"Data/Microarray/analysis_output/WGCNA/modules.rds"


SOMALOGIC_STABILITY_IN_PATH <- snakemake@input[["soma_stability"]]
MARRAY_STABILITY_IN_PATH <- snakemake@input[["array_stability"]]

SOMALOGIC_ESET_IN_PATH <- snakemake@input[["soma_eset"]]
MARRAY_ESET_IN_PATH <- snakemake@input[["array_eset"]]

SOMALOGIC_OUT_PATH <- snakemake@output[["soma"]]#"Paper_1_Figures/Figure_1_Tables/somalogic_module_members.csv"
MICROARRAY_OUT_PATH <- snakemake@output[["array"]]#"Paper_1_Figures/Figure_1_Tables/microarray_module_members.csv"

SOMALOGIC_FEATCOUNT_OUT_PATH <- snakemake@output[["soma_feat_counts"]]#"Paper_1_Figures/Figure_1_Tables/somalogic_module_members.csv"
MICROARRAY_FEATCOUNT_OUT_PATH <- snakemake@output[["array_feat_counts"]]#"Paper_1_Figures/Figure_1_Tables/microarray_module_members.csv"
source("scripts/util/Plotting/tbnk_featurename_replace.R")

soma_mods <- readRDS(SOMALOGIC_MODULE_IN_PATH)

array_mods <- readRDS(MICROARRAY_MODULE_IN_PATH)

soma_stability <- readRDS(SOMALOGIC_STABILITY_IN_PATH)
array_stability <- readRDS(MARRAY_STABILITY_IN_PATH)

soma_stab_feat <- rownames(soma_stability)
array_stab_feat <- rownames(array_stability)

soma_eset <- readRDS(SOMALOGIC_ESET_IN_PATH)
array_eset <- readRDS(MARRAY_ESET_IN_PATH)

soma_featdata <- featureData(soma_eset)@data
array_featdata <- featureData(array_eset)@data

array_featdata <- array_featdata %>%
        rownames_to_column(var = "feature") %>%
        .[, !grepl("UniGene", colnames(.))] %>%
        .[, !grepl("GO", colnames(.))] %>%
        .[, !grepl("Chromosome", colnames(.))] %>%
        select(-c(Nucleotide.Title, Platform_CLONEID, Platform_ORF))

soma_featdata <- soma_featdata %>%
        rownames_to_column(var = "feature") %>%
        select(-c(Units, Type, Dilution))


soma_mods_dat <-
        soma_mods %>%
        enframe(value = "module_color", name = "feature") %>%
        mutate(module_name = replace_mod_names_single_type(module_color, sheet = "PM")) %>%
        mutate(tmp = gsub("grey", "PM9999", module_name)) %>%
        left_join(soma_featdata) %>%
        mutate(stable = feature %in% soma_stab_feat) %>%
        arrange(tmp) %>%
        select(-tmp)

soma_mods_dat %>%
        write_csv(path = SOMALOGIC_OUT_PATH)

array_mods_dat <- array_mods %>%
        enframe(value = "module_color", name = "feature") %>%
        mutate(module_name = replace_mod_names_single_type(module_color, sheet = "TM")) %>%
        left_join(array_featdata) %>%
        mutate(stable = feature %in% array_stab_feat) %>%
        arrange(module_name)

array_mods_dat %>%
        write_csv(path = MICROARRAY_OUT_PATH)


array_feat_counts <- array_mods_dat %>%
        group_by(module_name, module_color, stable) %>%
        summarise(n = n()) %>%
        ungroup %>%
        spread(key = stable, value = n) %>%
        rename(n_stable = "TRUE", n_unstable = "FALSE") %>%
        mutate(n_total = n_stable + n_unstable)

soma_feat_counts <- soma_mods_dat %>%
        group_by(module_name, module_color, stable) %>%
        summarise(n = n()) %>%
        ungroup %>%
        spread(key = stable, value = n) %>%
        rename(n_stable = "TRUE", n_unstable = "FALSE") %>%
        mutate(n_total = n_stable + n_unstable)


write_csv(array_feat_counts, MICROARRAY_FEATCOUNT_OUT_PATH)
write_csv(soma_feat_counts, SOMALOGIC_FEATCOUNT_OUT_PATH)
