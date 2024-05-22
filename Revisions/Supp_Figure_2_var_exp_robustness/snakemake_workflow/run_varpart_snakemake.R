library(variancePartition)
library(BiocParallel)
library(SummarizedExperiment)
library(edgeR)
library(readr)
#library(doParallel)

#cl <- makeCluster(4)
#registerDoParallel(cl)
register(SerialParam())
registered()

# Input Paths -----------------------------------------------------
ESET.IN.PATH <- snakemake@input[[1]]

# Output Paths ----------------------------------------------------
VARPART.OUT.PATH <- snakemake@output[[1]]
META.OUT.PATH <- snakemake@output[[2]]
SUBJ.REP.OUT.PATH <- snakemake@output[[3]]
SUBJ.NOREP.OUT.PATH <- snakemake@output[[4]]

# Read in data ----------------------------------------------------
eset <- readRDS(ESET.IN.PATH)


boot <- snakemake@params[["boot"]]
set.seed(boot)

subj_w_repeats <- unique(eset$patient_id[duplicated(eset$patient_id)])
subj_wo_repeats <- setdiff(unique(eset$patient_id), subj_w_repeats)

n_subj_rep_keep <- round(length(subj_w_repeats) * .8)
n_subj_norep_keep <- round(length(subj_wo_repeats) * .8)
keep_subj_rep <- sample(as.character(subj_w_repeats), size = n_subj_rep_keep)
keep_subj_norep <- sample(as.character(subj_wo_repeats), size = n_subj_norep_keep)

keep_subj_combined <- c(keep_subj_rep, keep_subj_norep)

eset <- eset[, eset$patient_id %in% keep_subj_combined]

meta <- pData(eset)
saveRDS(list(keep_subj_norep, keep_subj_rep, meta), "soma_meta_test_list.rds")
mat <- exprs(eset)

writeLines(as.character(keep_subj_rep), SUBJ.REP.OUT.PATH)
writeLines(as.character(keep_subj_norep), SUBJ.NOREP.OUT.PATH)

# Define formula --------------------------------------------------
form <-
  ~ (1|patient_id) 

# Run variancePartition analysis ----------------------------------
varPart <- fitExtractVarPartModel(mat, form, meta)

saveRDS(varPart, VARPART.OUT.PATH)
write_tsv(meta, META.OUT.PATH)
