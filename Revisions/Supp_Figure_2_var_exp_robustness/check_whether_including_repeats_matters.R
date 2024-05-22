library(variancePartition)
library(BiocParallel)
library(dplyr)
register(SerialParam())
registered()

tbnk <- readRDS("../../2021_07_08/Pipeline_out/Data/TBNK/processed/tbnk_eset_training.rds")


keepfeat <- c("cd4_cd3", "cd8_cd3", "cd19", "nk_cells", 
              "neutrophil_percent", "monocytes_percent", "eosinophil_percent",
              "basophil_percent")

tbnk <- tbnk[keepfeat, ]
tbnk <- tbnk[, complete.cases(t(exprs(tbnk)))]


tbnk_mat <- exprs(tbnk)
meta <- pData(tbnk)
meta$patient_id <- factor(meta$patient_id)

form <-
  ~ (1|patient_id) 

vp1 <- fitExtractVarPartModel(tbnk_mat, form, meta)

dup_pats <- meta$patient_id %>% .[duplicated(.)] %>% unique() %>% as.character()

keep_samp <- meta$patient_id %in% dup_pats


vp2 <- fitExtractVarPartModel(tbnk_mat[, keep_samp], form, meta[keep_samp, ])

vp1
vp2
