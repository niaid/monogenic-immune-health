#https://stackoverflow.com/questions/59645426/error-when-exporting-r-data-frame-using-openxlsx-error-in-zipr
#I was getting the same error.
#
#I think the problem is with dependencies of openxlsx. There is a "zipR" package that might be picked up when you install openxlsx, while the actual dependency is zip package:
#
#https://cran.r-project.org/web/packages/zip/index.html
#https://cran.r-project.org/web/packages/zipR/zipR.pdf
#I installed "zip" along with openxlsx and I don't get the error anymore

library(tidyverse)
library(data.table)
library(zip)
library(openxlsx)

#setwd("/hpcdata/sg/sg_data/PROJECTS/Monogenic_Project")

source("scripts/util/paper/abbrev_cond.R")

if(!exists("snakemake")){
  setwd("../..")
  source("scripts/util/paper/parse_snakemake.R")
  parse_snakemake("compile_supp_tables_excel")
}


paths_list <- list(
  s2_TM_members=snakemake@input[["tm_members"]],#"Paper_1_Figures/Figure_1_Tables/microarray_module_members.csv",
  s3_PM_members=snakemake@input[["pm_members"]],#"Paper_1_Figures/Figure_1_Tables/somalogic_module_members.csv",
  s4_TM_enrich=snakemake@input[["tm_enrich"]],#"Paper_1_Figures/Figure_1_Tables/microarray_module_gene_set_enrichments_table.txt",
  s5_PM_enrich_geneset=snakemake@input[["pm_enrich_genesets"]],#"Paper_1_Figures/Figure_1_Tables/somalogic_module_gene_set_enrichments_table.txt",
  s6_PM_enrich_tissue=snakemake@input[["pm_enrich_tissue"]],#"Paper_1_Figures/Figure_1_Tables/somalogic_module_tissues_set_enrichments_table.txt",
  s7_tbnk_DE=snakemake@input[["tbnk_de"]],#"Paper_1_Figures/Figure_2_Tables/cbc_and_tbnks_DE_results.txt",
  s8_PM_DE=snakemake@input[["pm_de"]],#"Paper_1_Figures/Figure_2_Tables/protein_modules_DE_results.txt",
  S9_TM_DE=snakemake@input[["tm_de"]],#"Paper_1_Figures/Figure_2_Tables/gene_modules_DE_results.txt",
  s10_P_feat_DE=snakemake@input[["p_feat_de"]],#"Paper_1_Figures/Figure_2_Tables/protein_features_DE_results.txt",
  s11_T_feat_DE=snakemake@input[["t_feat_de"]],#"Paper_1_Figures/Figure_2_Tables/gene_features_DE_results.txt",
  s12_Jive_PCs=snakemake@input[["jpcs"]],#"Paper_1_Figures/Figure_3_Tables/jive_pcs.csv",
  s13_JIVE_PC_cor_feat=snakemake@input[["jive_pc_feat_cor"]],
  s14_Jive_PC_enrich=snakemake@input[["jpc_enrich"]],#"Paper_1_Figures/Figure_3_Tables/jive_pc_enrichment.csv",
  s15_JIVE_PC_cor_mod_tbnk=snakemake@input[["jpc_cor_tbnk"]],
  s16_IHM_feat_gvi=snakemake@input[["ihm_feat_gvi"]],#"Paper_1_Figures/Figure_4_Tables/healthy_feature_gvi_table.txt",
  s17_IHM_scores=snakemake@input[["ihm_scores"]],#"Paper_1_Figures/Figure_4_Tables/hi_results_full_mod.csv",
  s18_meta_analysis_subjects=snakemake@input[["meta_analysis_n_subj"]],
  s19_sig_genes=snakemake@input[["sig_genes"]],#"Paper_1_Figures/Figure_4_Tables/surrogate_sig_genes.csv",
  s20_meta_analysis=snakemake@input[["meta_analysis"]],#"Paper_1_Figures/Figure_4_Tables/figure_4_meta_analysis_table.txt",
  s21_study_eff_size=snakemake@input[["study_eff_size"]],#"Paper_1_Figures/Figure_4_Tables/figure_4_meta_analysis_effsize_table.txt",
  s22_IHM_sig_prot=snakemake@input[["ihm_sig_prot"]],#"Paper_1_Figures/Figure_5_Tables/proteomic_surrogate_ihm.csv"
  s23_ihm_age_cxcl9_lm=snakemake@input[["ihm_age_cxcl9"]]
)

EXCEL.OUT.PATH <- snakemake@output[[1]]#"Paper_1_Figures/supp_tables.xlsx"

dat_list <- lapply(paths_list, fread, data.table = FALSE, nThread = 1)

anno_list <- list(
  s2_TM_members="Supplementary Table 2 : Transcriptomic Module (TM) gene membership",
  s3_PM_members="Supplementary Table 3 : Proteomic Module (PM) gene membership",
  s4_TM_enrich="Supplementary Table 4 : Transcriptomic Module (TM) gene set enrichment",
  s5_PM_enrich_geneset="Supplementary Table 5 : Proteomic Module (PM) gene set enrichment",
  s6_PM_enrich_tissue="Supplementary Table 6 : Proteomic Module (PM) tissue enrichment",
  s7_tbnk_DE="Supplementary Table 7 : CBC + TBNK feature differential abundance",
  s8_PM_DE="Supplementary Table 8 : Proteomic Module (PM) differential expression",
  S9_TM_DE="Supplementary Table 9 : Transcriptomic Module (TM) differential expression",
  s10_P_feat_DE="Supplementary Table 10 : Proteomic feature differential expression",
  s11_T_feat_DE="Supplementary Table 11 : Transcriptomic feature differential expression",
  s12_Jive_PCs="Supplementary Table 12 : JIVE Principal Component scores for each subject",
  s13_JIVE_PC_cor_feat="Supplementary Table 13 : JIVE Principal Component scores correlation with transcriptomic and proteomic features",
  s14_Jive_PC_enrich="Supplementary Table 14 : JIVE Principal Component Enrichment",
  s15_JIVE_PC_cor_mod_tbnk="Supplementary Table 15 : JIVE Principal Component scores correlation with modules and cell frequencies",
  s16_IHM_feat_gvi="Supplementary Table 16 : Classifier Global Variable Importance (GVI)",
  s17_IHM_scores="Supplementary Table 17 : Immune Health Metric (IHM) scores",
  s18_meta_analysis_subjects="Supplementary Table 18 : Studies included in meta-analysis and number of subjects",
  s19_sig_genes="Supplementary Table 19 : Transcriptomic surrogate signature genes used in meta-analysis",
  s20_meta_analysis="Supplementary Table 20 : Meta-analysis summary",
  s21_study_eff_size="Supplementary Table 21 : Meta-analysis within study effect sizes for the IHM surrogate signature",
  s22_IHM_sig_prot="Supplementary Table 22 : Proteomic Immune Health Metric (IHM) surrogate proteins",
  s23_ihm_age_cxcl9_lm="Supplementary Table 23 : Linear model results: IHM ~ age + cxcl9. In Monogenic and Baltimore data."
)

anno_list2 <- list(
  s2_TM_members="Additional Details: Modules were created with WGCNA R package. Both stable and unstable features were included. A feature's stability (variance explained by subject > .5) is shown. The number of features in each module are also shown.",
  s3_PM_members="Additional Details: Modules were created with WGCNA R package. Both stable and unstable features were included. A feature's stability (variance explained by subject > .5) is shown. The number of features in each module are also shown.",
  s4_TM_enrich="Additional Details: Enrichment determined with Fisher's Exact Test. P values adjusted with Benjamini-Hochberg procedure.",
  s5_PM_enrich_geneset="Additional Details: Enrichment determined with Fisher's Exact Test. P values adjusted with Benjamini-Hochberg procedure.",
  s6_PM_enrich_tissue="Additional Details: '25 Data from the Human Protein Atlas in tab-separated format', proteinatlas.tsv, was downloaded from https://www.proteinatlas.org/about/download. Genes were grouped into 3 categories per tissue based on characterization by Human Protein Atlas. strict = 'Tissue enriched', medium = 'Tissue enriched' + 'Tissue enhanced', general = 'Tissue enhanced + 'Tissue enriched' + 'Group enriched'. See source column for description of categories were included for that particular set.",
  s7_tbnk_DE="Additional Details: Linear models fit with Limma R package, comparing each condition to healthy controls while controlling for age and sex.",
  s8_PM_DE="Additional Details: Linear models fit with Limma R package, comparing each condition to healthy controls while controlling for age and sex.",
  S9_TM_DE="Additional Details: Linear models fit with Limma R package, comparing each condition to healthy controls while controlling for age and sex.",
  s10_P_feat_DE="Additional Details: Linear models fit with Limma R package, comparing each condition to healthy controls while controlling for age and sex.",
  s11_T_feat_DE="Additional Details: Linear models fit with Limma R package, comparing each condition to healthy controls while controlling for age and sex.",
  s12_Jive_PCs="Additional Details: Data were averaged for each subject and the stable transcriptomic and serum protein features were selected and used to to compute JIVE PC scores.",
  s13_JIVE_PC_cor_feat="Additional Details: The JIVE PC Scores for every subject were tested for correlation with all modules and cell population frequenceis",
  s14_Jive_PC_enrich="Additional Details: The correlation of each gene in the whole blood transcriptome data was computed. Gene set enrichment was then performed with the CameraPR function from the Limma R package",
  s15_JIVE_PC_cor_mod_tbnk="Additional Details: The JIVE PC Scores for every subject were tested for correlation with all modules and cell population frequenceis.",
  s16_IHM_feat_gvi="Additional Details: The Global Variable Importance is a measure of how useful a particular feature was to the classifier. P values were determined through a permutation test",
  s17_IHM_scores="Additional Details: IHM scores are the leave one out cross-validation scores predicting Healthy vs. Disease for each subject using the TMs, PMs, cell frequencies and grey module proteins. A higher score indicates a that this subject is more similar to the healthy subjects according to the classifier",
  s18_meta_analysis_subjects="Additional Details: Comparison Group Pairs described in Lau et al. F1000Research 5 (2016) were combined for each study and further curated as described in 'notes' column.",
  s19_sig_genes="Additional Details: Surrogate transcriptomic signatures of predictive features from IHM classifier were derived by searching for transcriptomic features highly correlated with the feature of interest (e.g. PM or cell population frequency)",
  s20_meta_analysis="Additional Details: Meta-analysis of transcriptomic surrogate signatures in autoimmunity datasets was performed with MetaIntegrator R package.",
  s21_study_eff_size="Additional Details: Meta-analysis of transcriptomic surrogate signature of IHM in autoimmunity datasets was performed with MetaIntegrator R package.",
  s22_IHM_sig_prot="Additional Details: Surrogate protein signature of IHM was derived by searching for protein features correlated with IHM",
  s23_ihm_age_cxcl9_lm="Additional Details: For the monogenic data, the IHM from the classifier was used directly in the linear model. For the Baltimore Aging cohort, the IHM proteomic surrogate was used."
)


wb <- createWorkbook()
## Add a worksheet
for(nm in names(dat_list)){
  dat <- dat_list[[nm]]
  if("condition" %in% colnames(dat)){
    dat$condition <- abbrev_cond(dat$condition)
  }
  if("adj.P.Val" %in% colnames(dat)){
    dat <- dat %>% rename(Adjusted.Pvalue = adj.P.Val)
  }
  if("p.adj" %in% colnames(dat)){
    dat <- dat %>% rename(Adjusted.Pvalue = p.adj)
  }
  if("AdjP" %in% colnames(dat)){
    dat <- dat %>% rename(Adjusted.Pvalue = AdjP)
  }
  if("p" %in% colnames(dat)){
    dat <- dat %>% rename(P.Value = p)
  }
  if("Pval" %in% colnames(dat)){
    dat <- dat %>% rename(P.Value = Pval)
  }
  colnames(dat) <- gsub("FDR", "AdjustedPVal", colnames(dat))

  addWorksheet(wb, nm)
  writeData(wb, nm, anno_list[[nm]])
  writeData(wb, nm, anno_list2[[nm]], startRow  = 3)
  writeData(wb, nm, dat, startRow = 5)
}



#Add featcounts for the modules for stable vs unstable
TM_FEAT_COUNT_IN_PATH <- snakemake@input[["tm_feat_counts"]]
PM_FEAT_COUNT_IN_PATH <- snakemake@input[["pm_feat_counts"]]
tm_feat_count <- fread(TM_FEAT_COUNT_IN_PATH, data.table = FALSE, nThread = 1)
pm_feat_count <- fread(PM_FEAT_COUNT_IN_PATH, data.table = FALSE, nThread = 1)
writeData(wb, 1, tm_feat_count, startRow = 5, startCol = ncol(dat_list[[1]]) + 3)
writeData(wb, 2, pm_feat_count, startRow = 5, startCol = ncol(dat_list[[2]]) + 3)

saveWorkbook(wb, file = EXCEL.OUT.PATH, overwrite = TRUE)

