library(tidyverse)
library(Biobase)
library(ggpubr)

#if(exists("snakemake")){
AGING.ESET.IN.PATH = snakemake@input[["aging_eset"]]#"Reference/ferrucci/processed/aging_eset.RDS"
PROTEOMIC.SIGNATURE.IN.PATH = snakemake@input[["proteomic_surrogate"]]#'Classification/proteomic_surrogates/healthy.index.plasma.surrogate.somaId.RDS'

PLOT.OUT.PATH = snakemake@output[[1]]#"Paper_1_Figures/Figure_5/il6_cor_immune_health_metric_in_ferrucci_using_surrogate.pdf"


#} else {
#  AGING.ESET.IN.PATH = "Reference/ferrucci/processed/aging_eset.RDS"
#  PROTEOMIC.SIGNATURE.IN.PATH = 'Classification/proteomic_surrogates/healthy.index.plasma.surrogate.somaId.RDS'
#
#
#  PLOT.OUT.PATH = "Paper_1_Figures/Figure_5/il6_cor_immune_health_metric_in_ferrucci_using_surrogate.pdf"
#
#  setwd("/hpcdata/sg/sg_data/PROJECTS/Monogenic_Project")
#}

source('scripts/util/Signatures/get_signature_scores.R')
# Figure 5a -- Baltimore Healthy Aging Study, age versus HI surrogate signature scatterplot with regression line

## Load the Baltimore Healthy Aging Study eset and the somalogic Healthy Index plasma surrogate signature
eset = readRDS(AGING.ESET.IN.PATH)
healthy.index.surrogate.signature = readRDS(PROTEOMIC.SIGNATURE.IN.PATH)

## Subset the eset to just the samples (not QC, Calibrators, or buffers)
eset = eset[,eset$SampleType == 'Sample']

## Extract the scores from the signature
X = t(exprs(eset))
scores = util.get_signature_score(X, healthy.index.surrogate.signature)

soma_mat <- exprs(eset)
feat_names <- featureData(eset)@data
grep("il-6", feat_names$Target, ignore.case = T, value = T)

il6_id <- feat_names %>%
        filter(Target == "IL-6") %>% 
        pull(SomaId)

il6_scores <- soma_mat[il6_id, ]

## Put the healthy index surrogate scores into a data frame
df = data.frame(`IL-6`= il6_scores, `Immune Health Metric`= scores, check.names = FALSE)

p <- ggplot(df, aes(x = `IL-6`, y = `Immune Health Metric`)) +
        geom_point() +
        stat_cor(method = "spearman") +
        theme_bw() +
        ggtitle("Ferrucci data using IHM surrogate")

ggsave(plot = p, filename = PLOT.OUT.PATH, height = 3, width = 3)


