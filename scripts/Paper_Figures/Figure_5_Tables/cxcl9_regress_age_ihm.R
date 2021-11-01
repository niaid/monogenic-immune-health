library(tidyverse)
library(Biobase)
library(BiocGenerics)

if(!exists("snakemake")){
  setwd("../../..")
  source("scripts/util/paper/parse_snakemake.R")
  parse_snakemake(rule = "figure_5_tables_cxcl9_ihm_age_regression")
}


AGING.ESET.IN.PATH = snakemake@input[[1]]
PROTEOMIC.SIGNATURE.IN.PATH = snakemake@input[[2]]

PLOT.OUT.PATH = snakemake@output[[1]]#"Paper_1_Figures/Figure_5/cxcl9_cor_immune_health_metric_in_ferrucci_using_surrogate.pdf"


source('scripts/util/Signatures/get_signature_scores.R')
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
#grep("il-6", feat_names$Target, ignore.case = T, value = T)
grep("MIG", feat_names$Target, ignore.case = T, value = T)

cxcl9_id <- feat_names %>%
        filter(Target == "MIG") %>% 
        pull(SomaId)

cxcl9_scores <- soma_mat[cxcl9_id, ]

## Put the healthy index surrogate scores into a data frame
df = data.frame(`cxcl9`= cxcl9_scores, `ihm`= scores,
                age = eset$Age,
                check.names = FALSE)

bmore_mod <- lm(ihm ~ age + cxcl9, data = df)
bmore_summ <- summary(bmore_mod)
bmore_dat <- as.data.frame(bmore_summ$coefficients) %>%
        rownames_to_column(var = "term")



#---- monogenic
hi.in.path = snakemake@input[[3]]
meta.in.path = snakemake@input[[4]]
soma.in.path = snakemake@input[[5]]

#fig.out.path = snakemake@output[[1]]#"paper_1_figures/figure_5/cxcl9_cor_immune_health_metric.pdf"

hi_dat <- readRDS(hi.in.path)
meta <- readRDS(meta.in.path)
soma <- readRDS(soma.in.path)


hi_dat <- hi_dat %>%
        rownames_to_column(var = "patient_id") %>%
        mutate(`ihm` = all.modules.plus.grey.with.tbnks) %>%
        select(patient_id, `ihm`)

soma_mat <- exprs(soma)
grep("9", rownames(soma_mat), ignore.case = t, value = t)
grep("mig", rownames(soma_mat), ignore.case = t, value = t)

feat_dat <- featureData(soma)@data

feat_dat %>% filter(Target == "MIG")

cxcl9_dat <- pData(soma) %>%
        select(patient_id, condition) %>%
        mutate(`cxcl9` = soma_mat["MIG", ]) %>%
        mutate(age = soma$Age)

df <- left_join(hi_dat, cxcl9_dat)

## put the healthy index surrogate scores into a data frame
#df = data.frame(`cxcl9`= , `ihm`= scores,
#                age = eset$age,
#                check.names = false)

mono_mod <- lm(ihm ~ age + cxcl9, data = df)
mono_summ <- summary(mono_mod)
mono_dat <- as.data.frame(mono_summ$coefficients) %>%
        rownames_to_column(var = "term")


combined_dat <- bind_rows(list(Monogenic = mono_dat, Baltimore = bmore_dat), .id = "Study")

TAB.OUT.PATH <- snakemake@output[[1]]

write_csv(combined_dat, TAB.OUT.PATH)
