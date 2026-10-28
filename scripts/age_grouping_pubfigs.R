# scripts/age_grouping_pubfigs.R
# Produce publication-quality limma heatmap (PNG 300dpi and PDF)

library(pheatmap)
library(RColorBrewer)
library(viridis)

DATA_FILE <- "data/expr_sample_table.csv"  # set to your real SomaScan/Olink CSV
PROTEIN_LIST_FILE <- "data/top50_proteins.txt" # optional: one protein per line
LIMMA_RES <- "outputs/limma_results_trend.csv"
OUT_DIR <- "outputs"
TOPN <- 50
PROT_PREFIX <- "p"  # set to '' to auto-detect numeric protein columns (recommended for SomaScan/Olink)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("Loading data...\n")
df <- read.csv(DATA_FILE, stringsAsFactors = FALSE)
res <- read.csv(LIMMA_RES, stringsAsFactors = FALSE)

# proteins present in data
if(PROT_PREFIX == ""){
  # auto-detect numeric columns excluding common metadata
  meta_cols <- c('age','sex','batch','CRP','SampleID','age_q','age_grp','age_ord')
  prot_cols <- names(df)[sapply(df, is.numeric) & !(names(df) %in% meta_cols)]
} else {
  prot_cols <- grep(paste0('^', PROT_PREFIX), names(df), value = TRUE)
}
if(length(prot_cols) == 0) stop('No protein columns found')

# choose top N by adjusted p-value or use provided protein list
if(file.exists(PROTEIN_LIST_FILE)){
  listed <- readLines(PROTEIN_LIST_FILE)
  chosen <- intersect(trimws(listed), prot_cols)
  if(length(chosen) == 0) stop('Protein list provided but none found in data')
} else {
  if(!('adj.P.Val' %in% colnames(res))){
    # fallback common column name
    if('adj.P.Val' %in% colnames(res)){} else stop('limma results missing adj.P.Val')
  }
  res_ordered <- res[order(res$adj.P.Val), ]
  chosen <- head(res_ordered$protein, n = TOPN)
  chosen <- intersect(chosen, prot_cols)
}

if(length(chosen) == 0) stop('No chosen proteins found in data')

# expression matrix (rows=proteins, cols=samples)
expr_mat <- as.matrix(t(df[, prot_cols, drop = FALSE]))
colnames(expr_mat) <- if(!is.null(df$SampleID)) df$SampleID else paste0('S', seq_len(ncol(expr_mat)))

sub_expr <- expr_mat[chosen, , drop = FALSE]
# z-score per protein (row)
scaled <- t(scale(t(sub_expr)))

# annotations --- compute quartiles/bins if missing
if(!('age_q' %in% colnames(df))){
  q_breaks <- quantile(df$age, probs = seq(0,1,0.25), na.rm = TRUE)
  df$age_q <- cut(df$age, breaks = q_breaks, include.lowest = TRUE, labels = c('Q1','Q2','Q3','Q4'))
}
if(!('age_grp' %in% colnames(df))){
  mybreaks <- c(-Inf,50,65,80,Inf)
  df$age_grp <- cut(df$age, breaks = mybreaks, include.lowest = TRUE, labels = c('<50','50-64','65-79','80+'))
}
# ensure lengths match columns
ann_col <- data.frame(AgeQuartile = as.character(df$age_q), AgeBin = as.character(df$age_grp), Sex = as.character(df$sex), Batch = as.character(df$batch), CRP = df$CRP, stringsAsFactors = FALSE)
rownames(ann_col) <- colnames(scaled)

# color palette: use cividis (perceptually uniform, colorblind-friendly)
cols <- viridis::cividis(200)

# PDF (vector)
pdf_file <- file.path(OUT_DIR, 'age_diff_heatmap_limma_pub.pdf')
pdf(pdf_file, width = 10, height = 8)
p <- pheatmap(scaled,
              annotation_col = ann_col,
              show_rownames = TRUE,
              show_colnames = FALSE,
              clustering_method = 'complete',
              fontsize = 10,
              fontsize_row = 6,
              fontsize_col = 8,
              color = cols,
              border_color = NA)
invisible(dev.off())
cat('Wrote PDF:', pdf_file, '\n')

# PNG (high-res raster)
png_file <- file.path(OUT_DIR, 'age_diff_heatmap_limma_pub.png')
png(filename = png_file, width = 10, height = 8, units = 'in', res = 300)
# draw same heatmap
pheatmap(scaled,
         annotation_col = ann_col,
         show_rownames = TRUE,
         show_colnames = FALSE,
         clustering_method = 'complete',
         fontsize = 10,
         fontsize_row = 6,
         fontsize_col = 8,
         color = cols,
         border_color = NA)
invisible(dev.off())
cat('Wrote PNG:', png_file, '\n')

cat('Publication-quality heatmaps created in', OUT_DIR, '\n')
