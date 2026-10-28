# scripts/age_grouping_examples.R
# Usage: edit DATA_FILE to point to your CSV; then run with Rscript scripts/age_grouping_examples.R
# Expects CSV with rows = samples and columns: proteins named starting with 'p', and columns: age, sex, batch, CRP (optional), SampleID (optional)

library(limma)
library(pheatmap)
# Optional: library(sva)

set.seed(12345)

# === User params ===
DATA_FILE <- "data/expr_sample_table.csv"  # change to your path
OUT_DIR <- "outputs"
PROT_PREFIX <- "p"
TOPN <- 50
# ====================

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("Loading data:", DATA_FILE, "\n")
df <- read.csv(DATA_FILE, stringsAsFactors = FALSE)

# Age grouping: quartiles and custom
# quartiles
q_breaks <- quantile(df$age, probs = seq(0,1,0.25), na.rm = TRUE)
df$age_q <- cut(df$age, breaks = q_breaks, include.lowest = TRUE, labels = c("Q1","Q2","Q3","Q4"))
# custom bins example (edit as needed)
mybreaks <- c(-Inf, 50, 65, 80, Inf)
df$age_grp <- cut(df$age, breaks = mybreaks, include.lowest = TRUE, labels = c("<50","50-64","65-79","80+"))
# ordinal coding for trend
if(is.factor(df$age_q)) df$age_ord <- as.numeric(df$age_q) else df$age_ord <- as.numeric(factor(df$age_q))

# proteins
prot_cols <- grep(paste0('^', PROT_PREFIX), names(df), value = TRUE)
if(length(prot_cols) == 0) stop('No protein columns found with prefix ', PROT_PREFIX)

# expression matrix (rows=proteins, cols=samples)
expr_mat <- as.matrix(t(df[, prot_cols, drop = FALSE]))
colnames(expr_mat) <- if(!is.null(df$SampleID)) df$SampleID else paste0('S', seq_len(ncol(expr_mat)))

# design: use quartile factor (reference Q1) and covariates
# ensure factor ref
df$age_q <- relevel(factor(df$age_q), ref = "Q1")
design <- model.matrix(~ df$age_q + df$sex + df$batch + df$CRP)
colnames(design) <- make.names(colnames(design))

# limma fit
fit <- lmFit(expr_mat, design)
fit <- eBayes(fit)

# Find coef name for Q4 vs Q1 (example) or trend (age_ord)
# look for common age_q Q4 column name patterns
coef_name_q4 <- grep('age_qQ4', colnames(design), value = TRUE)
if(length(coef_name_q4) == 0) coef_name_q4 <- grep('Q4', colnames(design), value = TRUE)
if(length(coef_name_q4) == 0) coef_name_q4 <- NULL

# prefer trend test if present
coef_name_trend <- NULL
if('age_ord' %in% colnames(df)){
  design_trend <- model.matrix(~ df$age_ord + df$sex + df$batch + df$CRP)
  fit_trend <- lmFit(expr_mat, design_trend)
  fit_trend <- eBayes(fit_trend)
  # find numeric age_ord coefficient name in fit_trend
  coef_name_trend <- grep('age_ord', colnames(design_trend), value = TRUE)
  if(length(coef_name_trend)==0) coef_name_trend <- NULL
}

# If trend is available prefer it
if(!is.null(coef_name_trend)){
  res <- topTable(fit_trend, coef = coef_name_trend[1], number = Inf, adjust.method = 'BH')
  res$protein <- rownames(res)
  write.csv(res, file = file.path(OUT_DIR, 'limma_results_trend.csv'), row.names = FALSE)
  chosen <- head(rownames(res[order(res$adj.P.Val),]), n = TOPN)
} else {
  # attempt to find an age-related coefficient in the original design
  age_coefs <- grep('age_q|age', colnames(design), value = TRUE)
  age_coefs <- setdiff(age_coefs, '(Intercept)')
  if(length(age_coefs) > 0){
    res <- topTable(fit, coef = age_coefs[1], number = Inf, adjust.method = 'BH')
    res$protein <- rownames(res)
    write.csv(res, file = file.path(OUT_DIR, 'limma_results_age_group.csv'), row.names = FALSE)
    chosen <- head(rownames(res[order(res$adj.P.Val),]), n = TOPN)
  } else {
    stop('No age coefficient found in design matrix')
  }
}

# subset and z-score per protein
sub_expr <- expr_mat[chosen, , drop = FALSE]
scaled <- t(scale(t(sub_expr)))

ann_col <- data.frame(AgeQuartile = df$age_q, AgeBin = df$age_grp, Sex = df$sex, Batch = df$batch, CRP = df$CRP)
rownames(ann_col) <- colnames(scaled)

# heatmap
pheatmap(scaled, annotation_col = ann_col, show_rownames = TRUE, show_colnames = FALSE, clustering_method = 'complete', filename = file.path(OUT_DIR, 'age_diff_heatmap_limma.png'), width = 10, height = 8)

cat('R script finished. Results in', OUT_DIR, '\n')
