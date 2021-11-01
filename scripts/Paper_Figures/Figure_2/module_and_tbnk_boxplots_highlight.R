library(tidyverse)
library(Biobase)
library(ggpubr)


TBNK_IN_PATH <- snakemake@input[["tbnk_data"]]#"Data/TBNK/data_analysis_ready/tbnk_training_sample_level_eset.rds"
PROTEIN_IN_PATH <- snakemake@input[["soma_mod_scores"]]#"Data/Somalogic/analysis_output/wgcna_results/scores_sample_level.rds"
GENE_IN_PATH <- snakemake@input[["array_mod_scores"]]#"Data/Microarray/analysis_output/WGCNA/array_sample_scores.rds"


RESULTS.IN.PATHS = list(
  somalogic.modules = snakemake@input[["soma_mod_de"]],# 'Data/Somalogic/analysis_output/differential_expression/somalogic_modules_DE_results.rds',
  somalogic.features = snakemake@input[["soma_feat_de"]],# 'Data/Somalogic/analysis_output/differential_expression/somalogic_features_DE_results.rds',
  microarray.modules = snakemake@input[["array_mod_de"]],# 'Data/Microarray/analysis_output/differential_expression/microarray_modules_DE_results.rds',
  microarray.features = snakemake@input[["array_feat_de"]],# 'Data/Microarray/analysis_output/differential_expression/microarray_features_DE_results.rds',
  tbnks = snakemake@input[["tbnk_de"]]# 'Data/TBNK/analysis_output/differential_expression/tbnks_features_DE_results.rds'
)

FIG.OUT.PATH <- snakemake@output[[1]]#"Paper_1_Figures/Figure_2/boxplots_of_select_features_v2.pdf"
#setwd("../../..")
#TBNK_IN_PATH <- "Pipeline_out/Data/TBNK/data_analysis_ready/tbnk_training_sample_level_eset.rds"
#PROTEIN_IN_PATH <- "Pipeline_out/Data/Somalogic/analysis_output/wgcna_results/scores_sample_level.rds"
#GENE_IN_PATH <- "Pipeline_out/Data/Microarray/analysis_output/WGCNA/array_sample_scores.rds"
#RESULTS.IN.PATHS = list(
#  somalogic.modules = 'Pipeline_out/Data/Somalogic/analysis_output/differential_expression/somalogic_modules_DE_results.rds',
#  somalogic.features = 'Pipeline_out/Data/Somalogic/analysis_output/differential_expression/somalogic_features_DE_results.rds',
#  microarray.modules = 'Pipeline_out/Data/Microarray/analysis_output/differential_expression/microarray_modules_DE_results.rds',
#  microarray.features = 'Pipeline_out/Data/Microarray/analysis_output/differential_expression/microarray_features_DE_results.rds',
#  tbnks = 'Pipeline_out/Data/TBNK/analysis_output/differential_expression/tbnks_features_DE_results.rds'
#)
#FIG.OUT.PATH <- "Pipeline_out/Paper_1_Figures/Figure_2/boxplots_of_select_features_v2.pdf"

source("scripts/util/paper/abbrev_cond.R")

tbnk_eset <- readRDS(TBNK_IN_PATH)
protein_eset <- readRDS(PROTEIN_IN_PATH)
gene_eset <- readRDS(GENE_IN_PATH)

# Load data
#esets = lapply(ESETS.IN.PATHS, readRDS) # Load expression set data
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



results_q_dat <- results$all$versus.healthy$adj.P.Val %>%
        as.data.frame() %>%
        rownames_to_column(var = "feature") %>%
        gather(key = "condition", value = "q_val", -feature)


asterisk_vec <- function(p_vec){
  out_vec <- rep("", length(p_vec))
  out_vec <- replace(out_vec, p_vec < .05, "*")
  out_vec <- replace(out_vec, p_vec < .01, "**")
  out_vec <- replace(out_vec, p_vec < .001, "***")
  out_vec
}

nk_cell_signif <- results_q_dat %>%
        filter(feature =="tbnks.nk_cells_abs") %>%
        mutate(nk_cell_signif = asterisk_vec(q_val))

rdw_signif <- results_q_dat %>%
        filter(feature =="tbnks.rdw") %>%
        mutate(rdw_signif = asterisk_vec(q_val))

protein_magenta_signif <- results_q_dat %>%
        filter(feature =="protein.magenta") %>%
        mutate(protein_magenta_signif = asterisk_vec(q_val))

protein_purple_signif <- results_q_dat %>%
        filter(feature =="protein.purple") %>%
        mutate(protein_purple_signif = asterisk_vec(q_val))

gene_yellow_signif <- results_q_dat %>%
        filter(feature =="gene.yellow") %>%
        mutate(gene_yellow_signif = asterisk_vec(q_val))

gene_magenta_signif <- results_q_dat %>%
        filter(feature =="gene.magenta") %>%
        mutate(gene_magenta_signif = asterisk_vec(q_val))

signif_dat_list <- list(nk_cell_signif, rdw_signif, protein_magenta_signif, gene_yellow_signif, protein_purple_signif, gene_magenta_signif)
signif_dat_list <- lapply(signif_dat_list, function(dat){
  dat[, c("condition", grep("signif", colnames(dat), value = TRUE))]
})

signif_dat <- Reduce(full_join, signif_dat_list)


tbnk_meta <- tbnk_eset %>%
        pData() %>%
        select(patient_id, gender, condition) %>%
        mutate(nk_cells_abs = exprs(tbnk_eset)["nk_cells_abs", ]) %>%
        mutate(rdw = exprs(tbnk_eset)["rdw", ]) %>%
        left_join(signif_dat) %>%
        mutate(cond_group = group_cond(condition)) %>%
        mutate(condition = factor(condition))

protein_meta <- protein_eset %>%
        pData() %>%
        select(patient_id, gender, condition) %>%
        mutate(PM6 = exprs(protein_eset)["magenta", ]) %>%
        mutate(PM2 = exprs(protein_eset)["purple", ]) %>%
        left_join(signif_dat) %>%
        mutate(cond_group = group_cond(condition)) %>%
        mutate(condition = factor(condition))

gene_meta <- gene_eset %>%
        pData() %>%
        select(patient_id, gender, condition) %>%
        mutate(TM2= exprs(gene_eset)["yellow", ]) %>%
        mutate(TM6 = exprs(gene_eset)["magenta", ]) %>%
        left_join(signif_dat) %>%
        mutate(cond_group = group_cond(condition)) %>%
        mutate(condition = factor(condition))

levels(gene_meta$condition) <- abbrev_cond(levels(gene_meta$condition))
levels(protein_meta$condition) <- abbrev_cond(levels(protein_meta$condition))
levels(tbnk_meta$condition) <- abbrev_cond(levels(tbnk_meta$condition))


p1 <- ggplot(tbnk_meta, aes(y = nk_cells_abs, 
                            x = reorder(condition, nk_cells_abs, median))) +
        geom_boxplot(outlier.shape = NA, aes(fill = cond_group)) +
        ggbeeswarm::geom_beeswarm(alpha = .5, size = 1) +
        geom_text(aes(x = condition, y = max(nk_cells_abs) * 1.1, label = nk_cell_signif), size = 6, color = "red", position = position_nudge(x = -.5)) +
        ylim(min(tbnk_meta$nk_cells_abs),max(tbnk_meta$nk_cells_abs) * 1.2) +
        theme_bw() +
        ggtitle("NK cells (#)") +
        ylab("") + 
        coord_flip()+
        xlab("") + 
        facet_grid(condition == "Healthy" ~ 1, scales = "free_y", space = "free_y") +
        theme(legend.position = "none", 
              strip.background = element_blank(),
              strip.text.x = element_blank(), 
              strip.text.y = element_blank(),
              panel.spacing = unit(0, "lines"))

p2 <- ggplot(tbnk_meta, aes(y = rdw, x = reorder(condition, rdw , median)))+
                            #color = condition == "Healthy")) +
        #geom_boxplot(outlier.shape = NA) +
        geom_boxplot(outlier.shape = NA, aes(fill = cond_group)) +
        geom_text(aes(x = condition, y = max(rdw) * 1.1, label = rdw_signif), size = 6, color = "red", position = position_nudge(x = -.5)) +
        ylim(min(tbnk_meta$rdw),max(tbnk_meta$rdw) * 1.13) +
        #geom_jitter(width = 0, alpha = .5) +
        ggbeeswarm::geom_beeswarm(alpha = .5, size = 1) +
        theme_bw() +
        coord_flip()+
        ggtitle("RDW") +
        ylab("") + 
        xlab("") + 
        facet_grid(condition == "Healthy" ~ 1, scales = "free_y", space = "free_y") +
        theme(legend.position = "none", 
              strip.background = element_blank(),
              strip.text.x = element_blank(), 
              strip.text.y = element_blank(),
              panel.spacing = unit(0, "lines"))

p3 <- ggplot(gene_meta, aes(y = TM6, x = reorder(condition, TM6, median)))+
                            #color = condition == "Healthy")) +
        #geom_boxplot(outlier.shape = NA) +
        geom_boxplot(outlier.shape = NA, aes(fill = cond_group)) +
        geom_text(aes(x = condition, y = max(TM6) * 1.1, label = gene_magenta_signif), size = 6, color = "red", position = position_nudge(x = -.5)) +
        ylim(min(gene_meta$TM6),max(gene_meta$TM6) * 1.2) +
        #geom_jitter(width = 0, alpha = .5) +
        ggbeeswarm::geom_beeswarm(alpha = .5, size = 1) +
        xlab("TM6: CD8/NK cells") +
        coord_flip()+
        theme_bw() +
        facet_grid(condition == "Healthy" ~ 1, scales = "free_y", space = "free_y") +
        ggtitle("TM6: CD8/NK cells") +
        ylab("") + 
        xlab("") + 
        #scale_color_manual(values = c("black", "red")) + 
        theme(legend.position = "none", 
              strip.background = element_blank(),
              strip.text.x = element_blank(), 
              strip.text.y = element_blank(),
              panel.spacing = unit(0, "lines"))

p4 <- ggplot(gene_meta, aes(y = TM2, x = reorder(condition, TM2, median)))+
                            #color = condition == "Healthy")) +
        #geom_boxplot(outlier.shape = NA) +
        geom_boxplot(outlier.shape = NA, aes(fill = cond_group)) +
        geom_text(aes(x = condition, y = max(TM2) * 1.1, label = gene_yellow_signif), size = 6, color = "red", position = position_nudge(x = -.5)) +
        ylim(min(gene_meta$TM2),max(gene_meta$TM2) * 1.2) +
        #geom_jitter(width = 0, alpha = .5) +
        ggbeeswarm::geom_beeswarm(alpha = .5, size = 1) +
        xlab("TM2: heme/RBC score") +
        theme_bw() +
        coord_flip()+
        facet_grid(condition == "Healthy" ~ 1, scales = "free_y", space = "free_y") +
        ggtitle("TM2: heme/RBC") +
        xlab("") + 
        ylab("") + 
        #scale_color_manual(values = c("black", "red")) + 
        theme(legend.position = "none", 
              strip.background = element_blank(),
              strip.text.x = element_blank(), 
              strip.text.y = element_blank(),
              panel.spacing = unit(0, "lines"))

p5 <- ggplot(protein_meta, aes(y = PM6, x = reorder(condition, PM6, median))) +
                               #color = condition == "Healthy")) +
        #geom_boxplot(outlier.shape = NA) +
        geom_boxplot(outlier.shape = NA, aes(fill = cond_group)) +
        geom_text(aes(x = condition, y = max(PM6) * 1.1, label = protein_magenta_signif), size = 6, color = "red", position = position_nudge(x = -.5)) +
        ylim(min(protein_meta$PM6),max(protein_meta$PM6) * 1.22) +
        #geom_jitter(width = 0, alpha = .5) +
        ggbeeswarm::geom_beeswarm(alpha = .5, size = 1) +
        xlab("PM6: platelets score") +
        ggtitle("PM6: platelets") + 
        facet_grid(condition == "Healthy" ~ 1, scales = "free_y", space = "free_y") +
        coord_flip()+
        theme_bw() +
        ylab("") + 
        xlab("") + 
        #scale_color_manual(values = c("black", "red")) + 
        theme(legend.position = "none", 
              strip.background = element_blank(),
              strip.text.x = element_blank(), 
              strip.text.y = element_blank(),
              panel.spacing = unit(0, "lines"))

p6 <- ggplot(protein_meta, aes(y = PM2, x = reorder(condition, PM2, median))) +
                               #color = condition == "Healthy")) +
        #geom_boxplot(outlier.shape = NA) +
        geom_boxplot(outlier.shape = NA, aes(fill = cond_group)) +
        geom_text(aes(x = condition, y = max(PM2) * 1.1, label = protein_purple_signif), size = 6, color = "red", position = position_nudge(x = -.5)) +
        ylim(min(protein_meta$PM2),max(protein_meta$PM2) * 1.2) +
        #geom_jitter(width = 0, alpha = .5) +
        ggbeeswarm::geom_beeswarm(alpha = .5, size = 1) +
        xlab("PM2") +
        ggtitle("PM2") + 
        facet_grid(condition == "Healthy" ~ 1, scales = "free_y", space = "free_y") +
        coord_flip()+
        theme_bw() +
        ylab("") + 
        xlab("") + 
        #scale_color_manual(values = c("black", "red")) + 
        theme(legend.position = "none", 
              strip.background = element_blank(),
              strip.text.x = element_blank(), 
              strip.text.y = element_blank(),
              panel.spacing = unit(0, "lines"))


#file.remove(FIG.OUT.PATH)
#pdf(FIG.OUT.PATH, height = 9.6, width = 5.5)
#print(p1)
#print(p2)
#print(p3)
#print(p4)
#print(p5)
#print(p6)
#dev.off()
p <- cowplot::plot_grid(p1, p2, p3, p4, p5, p6, nrow = 3)
ggsave(plot =p, filename = FIG.OUT.PATH, height = 9.6, width = 5.7)
