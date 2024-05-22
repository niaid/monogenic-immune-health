
# Load libraries and source utilities function
library(tidyverse)
library(cowplot)

if(!exists("snakemake")){
  setwd("../../..")
  source("scripts/util/paper/parse_snakemake.R")
  parse_snakemake(rule = "figure_2_mod_feature_level_heatmaps")
}

# Set paths
RESULTS.IN.PATHS = list(
  somalogic.modules = snakemake@input[["soma_mod_de"]],# 'Data/Somalogic/analysis_output/differential_expression/somalogic_modules_DE_results.rds',
  somalogic.features = snakemake@input[["soma_feat_de"]],# 'Data/Somalogic/analysis_output/differential_expression/somalogic_features_DE_results.rds',
  microarray.modules = snakemake@input[["array_mod_de"]],# 'Data/Microarray/analysis_output/differential_expression/microarray_modules_DE_results.rds',
  microarray.features = snakemake@input[["array_feat_de"]],# 'Data/Microarray/analysis_output/differential_expression/microarray_features_DE_results.rds',
  tbnks = snakemake@input[["tbnk_de"]]# 'Data/TBNK/analysis_output/differential_expression/tbnks_features_DE_results.rds'
)

fdr_threshold <- .05

source("scripts/util/paper/abbrev_cond.R")


PROTEIN.INFO.IN.PATH = snakemake@input[["protein_info"]]

PROTEIN.MODULES.IN.PATH = snakemake@input[["soma_mods"]]#'Data/Somalogic/analysis_output/wgcna_results/modules.rds'
GENE.MODULES.IN.PATH = snakemake@input[["array_mods"]]#'Data/Microarray/analysis_output/WGCNA/modules.rds'
#PROTEIN.MODULES.IN.PATH = 'Data/Somalogic/analysis_output/wgcna_results/modules.rds'
#GENE.MODULES.IN.PATH = 'Data/Microarray/analysis_output/WGCNA/modules.rds'
COND.GROUPS.IN.PATH <- snakemake@input[["cond_groups"]]#"Reference/condition_groups.rds"

GENE.SETS.IN.PATH = snakemake@input[["genesets"]]#'Gene_sets/processed/combined_gene_sets.RDS'

PURPLE.MOD.OUT.PATH = snakemake@output[["purple_pm2"]]#"Paper_1_Figures/Figure_2/purple_protein_module_bubble_plot.pdf"
RED.TOP.ENRICH.OUT.PATH <- snakemake@output[["red_top_enrich"]]#"Paper_1_Figures/Figure_2/red_module_subcluster_enrich_top10.rds"
RED.MOD.CLUSTER.OUT.DIR <- snakemake@output[["red_subclus_dir"]]#"scratch/red_module_sub_clusters"
RED.DENDRO.OUT.PATH <- snakemake@output[["red_dendro"]]
RED.BUBBLE.OUT.PATH <- snakemake@output[["red_tm1"]]#"Paper_1_Figures/Figure_2/red_module_heatmap.pdf"

#
# Load data
#esets = lapply(ESETS.IN.PATHS, readRDS) # Load expression set data
results = lapply(RESULTS.IN.PATHS, readRDS) # Load results from limma DE testing

protein_info = read.csv(PROTEIN.INFO.IN.PATH, sep = "\t", stringsAsFactors = FALSE)

condition_groups <- readRDS(COND.GROUPS.IN.PATH)

gene.sets <- readRDS(GENE.SETS.IN.PATH)
#read in modules
gene_modules <- readRDS(GENE.MODULES.IN.PATH)
protein_modules <- readRDS(PROTEIN.MODULES.IN.PATH)

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
results_t_dat <- results$all$versus.healthy$t %>%
        as.data.frame() %>%
        rownames_to_column(var = "feature") %>%
        gather(key = "condition", value = "t_stat", -feature)

results_q_dat <- results$all$versus.healthy$adj.P.Val %>%
        as.data.frame() %>%
        rownames_to_column(var = "feature") %>%
        gather(key = "condition", value = "q_val", -feature)


joined_dat <- full_join(results_t_dat, results_q_dat) %>%
        mutate(q_val_filtered = replace(q_val, q_val > .2, NA)) %>%
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

#Plot the red module


prep_for_plotting <- function(res_list, keep_features = NULL, condition_groups){
  
  t_mat <- res_list$t
  coeff_mat <- res_list$effect.size
  q_mat <- res_list$adj.P.Val

  if(!is.null(keep_features)){
    coeff_mat <- coeff_mat[rownames(coeff_mat) %in% keep_features, ]
    t_mat <- t_mat[rownames(t_mat) %in% keep_features, ]
    q_mat <- q_mat[rownames(q_mat) %in% keep_features, ]
  }

  coeff_dat <- coeff_mat %>%
          as.data.frame() %>%
          rownames_to_column(var = "feature") %>%
          gather(key = "condition", value = "coeff", -feature)
  
  q_dat <- q_mat %>%
          as.data.frame() %>%
          rownames_to_column(var = "feature") %>%
          gather(key = "condition", value = "q_val", -feature)
  
  
  #joined <- full_join(t_dat, q_dat) #%>%
  #return(joined)
  joined <- full_join(coeff_dat, q_dat) %>%
          mutate(q_val_filtered = replace(q_val, q_val > .05, NA)) %>%
          left_join(condition_groups)

  
  hc_row <- t_mat %>% 
          dist() %>% hclust()
  joined$feature_ordered <- factor(joined$feature, levels = hc_row$labels[hc_row$order])
  
  hc_col <- t_mat %>% 
          t() %>%
          dist() %>% hclust()

  joined$condition_ordered <- factor(joined$condition, levels = hc_col$labels[hc_col$order])

  joined

}

red_genes <- names(gene_modules)[gene_modules == "red"]
red_dat <- prep_for_plotting(results$microarray.features$versus.healthy, 
                             keep_features = red_genes, condition_groups = condition_groups)

red_t_mat <- results$microarray.features$versus.healthy$t
red_t_mat <- red_t_mat[rownames(red_t_mat) %in% red_genes, ]

hc_red <- red_t_mat %>% 
        dist() %>%
        hclust()

pdf(RED.DENDRO.OUT.PATH)
plot(hc_red)
#abline(h = 12)
dev.off()

clusters_red <- cutree(hc_red, k = 3)

clusters_red_list <- lapply(unique(clusters_red), function(i){
  names(clusters_red)[clusters_red == i]                             
})

clusters_red_dat <- clusters_red_list %>% lapply(enframe, name = "name", value = "feature") %>% bind_rows(.id = "cluster")

red_dat <- left_join(red_dat, clusters_red_dat)

source("scripts/util/Enrichment/hyperGeo.R")

# Get set of all genes
universe = rownames(results$microarray.features$versus.healthy$t)


enrichments = lapply(clusters_red_list, function(hits) {
  multiHyperGeoTests(gene.sets, universe, hits, minGeneSetSize = 5, pAdjustMethod = "BH")
})


dir.create(RED.MOD.CLUSTER.OUT.DIR)
dir.create(file.path(RED.MOD.CLUSTER.OUT.DIR, "red_module_sub_clusters"))
for(i in seq_along(clusters_red_list)){
  path <- paste0(RED.MOD.CLUSTER.OUT.DIR, "/red_module_sub_clusters/genes_cluster_",i,".txt")
  writeLines(clusters_red_list[[i]], path)
  path2 <- paste0(RED.MOD.CLUSTER.OUT.DIR, "/red_module_sub_clusters/enrichments_cluster_",i,".csv")
  enrichments[[i]] %>% rownames_to_column(var = "geneset") %>%
  write_csv(path2)
}

top10_list <- lapply(enrichments, function(dat){
  dat %>% 
          rownames_to_column(var = "geneset") %>%
          arrange(across.Adjusted.Pvalue) %>% 
          head(10)
})
saveRDS(top10_list, RED.TOP.ENRICH.OUT.PATH)

show_features <- list(c1 = c(
                             "IFI16", "IFI27L2", "IFI35", "IFI44",
                             "IFI44L", "IFI6", "IFIH1", "IFIT1",  "IFIT2",
                             "IFIT3", "IFIT5", "IFITM1", "IFITM3", "IRF1", 
                             "IRF7", "IRF9", "ISG15"),
                      c2 = c("STAT1", "STAT2", "JAK2", "TRIM14", 
                             "TRIM5", "RIPK2", "POLB")
)

show_features <- unlist(show_features, use.names = FALSE)


ylabs <- replace(levels(red_dat$feature_ordered),
  !levels(red_dat$feature_ordered) %in% show_features,
  ""
)
names(ylabs) <- levels(red_dat$feature_ordered)

levels(red_dat$condition_ordered) <- abbrev_cond(levels(red_dat$condition_ordered))

pdf(RED.BUBBLE.OUT.PATH, height = 5, width =5)
p <- ggplot(red_dat %>% filter(!condition %in% c("CD14")), aes(x = condition_ordered, y = feature_ordered)) +
         geom_tile(aes(fill = coeff)) +
         #scale_color_viridis_c() +
         scale_fill_gradient2(low = "blue", high = "red") +
         facet_grid(cluster~cond_group, scales = "free", space = "free") +
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
               #strip.background.y = element_blank(), strip.text.y = element_blank(),
               axis.text.y = element_blank(), axis.ticks.y = element_blank())
print(p)

red_dat <- red_dat %>% mutate(showfeat = replace(feature_ordered, !feature_ordered %in% show_features, ""))

p <- ggplot(red_dat %>% filter(!condition %in% c("CD14")), aes(x = condition_ordered, y = feature_ordered)) +
         geom_tile(aes(fill = coeff)) +
         geom_text(aes(y = feature_ordered, x = 2, label = showfeat), size = 2) +
         #scale_color_viridis_c() +
         scale_fill_gradient2(low = "blue", high = "red") +
         facet_grid(cluster~cond_group, scales = "free", space = "free") +
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
               #strip.background.y = element_blank(), strip.text.y = element_blank(),
               axis.text.y = element_blank(), axis.ticks.y = element_blank())
print(p)
dev.off()



protein_info_sub <- protein_info %>%
        mutate(feature = make.names(Target)) %>%
        select(feature, Target)
tmp <- results$somalogic.features$versus.healthy

tmp <- lapply(tmp, function(x){
  rn <- protein_info_sub$Target[match(rownames(x), protein_info_sub$feature)]
  rownames(x) <- rn
  x
})

purple_genes <- names(protein_modules)[protein_modules == "purple"]
purple_genes <- protein_info_sub$Target[match(purple_genes, protein_info_sub$feature)]

purple_dat <- prep_for_plotting(tmp,
                                keep_features = purple_genes, 
                                condition_groups = condition_groups)

levels(purple_dat$feature_ordered)[levels(purple_dat$feature_ordered) == "MIG"] <-
        "CXCL9/MIG"

levels(purple_dat$condition_ordered) <- abbrev_cond(levels(purple_dat$condition_ordered))

pdf(PURPLE.MOD.OUT.PATH, height = 6, width = 6)
p <- ggplot(purple_dat %>% filter(!condition %in% c("TERT.TERC", "AI", "PID")), aes(x = condition_ordered, y = feature_ordered)) +
        #geom_point(aes(size = -log10(q_val), color = coeff)) +
        geom_point(aes(size = -log10(q_val), fill = coeff, color = q_val < fdr_threshold), shape = 21) +
        #scale_color_gradient2(low = "blue", high = "red") +
        scale_fill_gradient2(low = "blue", high = "red") +
        scale_color_manual(values = c("white", "black")) +
        theme_bw() +
        facet_grid(1~cond_group, scales = "free", space = "free") +
        ylab("") +
        theme(
              axis.text.y = element_text(size = 7.5),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7.5),
              strip.background = element_blank(), strip.text = element_blank())
print(p)

dev.off()

