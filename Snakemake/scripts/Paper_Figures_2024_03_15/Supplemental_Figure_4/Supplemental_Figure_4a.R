# Load libraries and source utilities function
library(Biobase)
library(tidyverse)
library(cowplot)

print(str(snakemake))
# Set paths
RESULTS.IN.PATHS = list(
  somalogic.modules = snakemake@input[["soma_mod_de"]],# 'Data/Somalogic/analysis_output/differential_expression/somalogic_modules_DE_results.rds',
  somalogic.features = snakemake@input[["soma_feat_de"]],# 'Data/Somalogic/analysis_output/differential_expression/somalogic_features_DE_results.rds',
  microarray.modules = snakemake@input[["array_mod_de"]],# 'Data/Microarray/analysis_output/differential_expression/microarray_modules_DE_results.rds',
  microarray.features = snakemake@input[["array_feat_de"]],# 'Data/Microarray/analysis_output/differential_expression/microarray_features_DE_results.rds',
  tbnks = snakemake@input[["tbnk_de"]]# 'Data/TBNK/analysis_output/differential_expression/tbnks_features_DE_results.rds'
)
#RESULTS.IN.PATHS = list(
#  somalogic.modules = 'Data/Somalogic/analysis_output/differential_expression/somalogic_modules_DE_results.rds',
#  somalogic.features = 'Data/Somalogic/analysis_output/differential_expression/somalogic_features_DE_results.rds',
#  microarray.modules = 'Data/Microarray/analysis_output/differential_expression/microarray_modules_DE_results.rds',
#  microarray.features = 'Data/Microarray/analysis_output/differential_expression/microarray_features_DE_results.rds',
#  tbnks = 'Data/TBNK/analysis_output/differential_expression/tbnks_features_DE_results.rds'
#)

fdr_threshold <- .05
source("scripts/util/paper/abbrev_cond.R")

COND.GROUPS.IN.PATH <- snakemake@input[["cond_groups"]]#"Reference/condition_groups.rds"

FIG.OUT.PATH <- snakemake@output[[1]]#"Paper_1_Figures/Figure_2/all_feature_bubble_plot_condition_separated.pdf"

condition_group_dat <- readRDS(COND.GROUPS.IN.PATH)

# Load data
results = lapply(RESULTS.IN.PATHS, readRDS) # Load results from limma DE testing


# Add a data type for the combination of protein modules, gene modules, and tbnks

## Get the names of the versus.healthy/versus.all options for the heatmaps
first.levels = names(results$somalogic.modules)
## Get the names of the statistics associated with the DE testing e.g. (t-stat, p value, etc.)
second.levels = names(results$somalogic.modules$versus.healthy)

## For both of the versus_healthy/versus_all options
results$all = lapply(first.levels, function(l1) {
  
  ## And for each statistic type associated with the DE testing results
  stats = lapply(second.levels, function(l2) {
    
    ## Get the matrix corresponding to the given statistic for the protein modules
    a1 = results$somalogic.modules[[l1]][[l2]]
    ## Name the protein modules with a 'protein.' prefix to distinguish them from the gene modules
    rownames(a1) = paste0('protein.', rownames(a1))
    ## Get the matrix corresponding to the given statistic for the gene modules
    a2 = results$microarray.modules[[l1]][[l2]]
    ## Name the gene modules with a 'gene.' prefix to distinguish them from the protein modules
    rownames(a2) = paste0('gene.', rownames(a2))
    ## Get the matrix corresponding to the given statistic for the tbnks
    a3 = results$tbnks[[l1]][[l2]]
    ## Name the tbnks with a 'tbnks.' prefix for consistency with the modules
    rownames(a3) = paste0('tbnks.', rownames(a3))
    
    ## Ensure the columns are in the same order
    a2 = a2[, colnames(a1)]
    a3 = a3[, colnames(a1)]
    
    ## Wrap these statistics together into a single matrix
    rbind(a1, a2, a3)
  })
  
  ## Name the statistics matrices with the corresponding statistic names (e.g. t-statistic, cohen's d, etc.)
  names(stats) = second.levels
  return(stats)
})

## Name the module-tbnk-combined results based on whether the versus.healthy option or versus.all option was used 
names(results$all) = first.levels

## Add the patient-condition relationships associated with all the esets in order to provide the number of samples for each condition
#sample.datas$all = unique(Reduce(rbind, sample.datas))

# Get the condition counts associated with each data type
#condition_counts = lapply(sample.datas, function(sample.data) {
#  counts = table(sample.data$condition)
#  counts = counts[names(counts) != "Healthy"]
#})



# make dataframe for first bubble_plot
joined_dat <- lapply(results$all, function(nested_list){
  lapply(names(nested_list), function(nm){
        mat <- nested_list[[nm]]
        dat <- as.data.frame(mat) %>%
        rownames_to_column(var = "feature") %>%
        gather(key = "condition", value = value, -feature)
        
         colnames(dat) <- gsub("value", nm, colnames(dat))
        dat
  }) %>% Reduce(f = full_join)
}) %>% bind_rows(.id = "comparison")


joined_dat <- joined_dat %>%
        rename(q_val = adj.P.Val, coeff = effect.size, t_stat = t) %>%
        mutate(q_val_filtered = replace(q_val, q_val > fdr_threshold, NA)) %>%
        mutate(feat_category = sapply(strsplit(feature, "\\."), `[[`, 1))


hc_row <-results$all$versus.healthy$t %>% 
        dist() %>% hclust()
joined_dat$feature <- factor(joined_dat$feature, levels = hc_row$labels[hc_row$order])

hc_col <-results$all$versus.healthy$t %>% 
        t() %>%
        dist() %>% hclust()
joined_dat$condition <- factor(joined_dat$condition, levels = hc_col$labels[hc_col$order])


#plotting everything all together
feat_types <- c("gene", "protein", "tbnks")
names(feat_types) <- feat_types
hc_row_list <- lapply(feat_types, function(feat_category){
  mat <- results$all$versus.healthy$t 
  mat[startsWith(rownames(mat), feat_category), ] %>% 
          dist() %>% hclust()
  
})

order_all_feat_within_datatype <- sapply(hc_row_list, function(hc){
  row_order <- hc$labels[hc$order]
}) %>% unlist()

joined_dat <- joined_dat %>%
        mutate(feature2 = factor(feature, levels = order_all_feat_within_datatype))

# separating out by feature type

lymph_terms <- c("cd19", "cd3", "nk_cell", "lymph")
innate_terms <- c("eos", "mono", "baso", "NLR", "wbc")
red_terms <- c("platelet", "mcv", "mch", "mchc", "rdw", "hemoglobin", "rbc")

joined_dat <- joined_dat %>%
        mutate(feat_group2 = feat_category) %>%
        mutate(feat_group2 = replace(feat_group2, 
                                     rowSums(sapply(lymph_terms, grepl, feature)) > 0,"Lymphocytes")
               ) %>%
        mutate(feat_group2 = replace(feat_group2, 
                                     rowSums(sapply(innate_terms, grepl, feature)) > 0,"Innate")
               ) %>%
        mutate(feat_group2 = replace(feat_group2, 
                                     rowSums(sapply(red_terms, grepl, feature)) > 0,"RBC & PLT")
               ) %>%
        mutate(feat_group2 = gsub("gene", "TM", feat_group2)) %>%
        mutate(feat_group2 = gsub("protein", "PM", feat_group2)) %>%
        mutate(feat_group2 = factor(feat_group2, levels = c("TM", "PM", "Innate", "Lymphocytes", "RBC & PLT")))

joined_dat <- joined_dat %>%
        mutate(cond_group = condition_group_dat$cond_group[match(condition, condition_group_dat$condition)])


#Add in the new module names that are not colors
joined_dat <- joined_dat %>%
        mutate(feature3 = feature2)

levels(joined_dat$feature3) <- 
        gsub("tbnks\\.", "", levels(joined_dat$feature3))

source("scripts/util/Plotting/tbnk_featurename_replace.R")

levels(joined_dat$feature3) <- replace_mod_names_both(levels(joined_dat$feature3))
levels(joined_dat$feature3) <- replace_tbnk_names(levels(joined_dat$feature3))

levels(joined_dat$condition) <- abbrev_cond(levels(joined_dat$condition))

pdf(FIG.OUT.PATH)
#Versus healthy- main figure
p <- ggplot(joined_dat %>% filter(comparison == "versus.healthy" & 
                                  !condition %in% c("TERT.TERC", "AI", "PID")), 
            aes(x = condition, y = feature3)) +
        geom_point(aes(size = -log10(q_val), fill = coeff, color = q_val < fdr_threshold), shape = 21) +
        facet_grid(feat_group2 ~ cond_group, scales = "free", space = "free") +
        scale_fill_gradient2(low = "blue", high = "red") +
        scale_color_manual(values = c("white", "black")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              strip.text.x = element_blank()) +
        ylab("") +
        ggtitle("Versus Healthy")
print(p)


#Versus all - this is supplement
p <- ggplot(joined_dat %>% filter(comparison == "versus.all" &
                                  !condition %in% c("TERT.TERC", "AI", "PID")), 
                                  aes(x = condition, y = feature3)) +
        geom_point(aes(size = -log10(q_val), fill = coeff, color = q_val < fdr_threshold), shape = 21) +
        facet_grid(feat_group2 ~ cond_group, scales = "free", space = "free") +
        theme_bw() +
        scale_fill_gradient2(low = "blue", high = "red") +
        scale_color_manual(values = c("white", "black")) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              strip.text.x = element_blank()) +
        ylab("") +
        ggtitle("Versus All")
print(p)
dev.off()

