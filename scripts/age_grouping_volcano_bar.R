# scripts/age_grouping_volcano_bar.R
# Create volcano plot and barplot of top effect sizes from limma results

library(ggplot2)
library(ggrepel)
library(viridis)

LIMMA_RES <- "outputs/limma_results_trend.csv"
OUT_DIR <- "outputs"
TOP_LABEL <- 50    # number of labels on volcano
TOPN_BAR <- 50     # number of bars in effect-size plot

res <- read.csv(LIMMA_RES, stringsAsFactors = FALSE)
if(!all(c('logFC','P.Value','adj.P.Val','protein','t') %in% colnames(res))) stop('limma results missing required columns')

# compute -log10 p
res$negLogP <- -log10(res$P.Value)
# significance threshold (adjusted)
res$signif <- res$adj.P.Val < 0.05
# fold-change threshold for highlighting (tunable)
fc_thr <- 0.8
res$highlight <- res$signif & abs(res$logFC) >= fc_thr

# continuous color by -log10 p for visual depth
point_cols <- viridis(100, option = "D")

# volcano plot
volcano_pdf <- file.path(OUT_DIR, 'limma_volcano_pub.png')
volcano_png <- file.path(OUT_DIR, 'limma_volcano_pub.png')

p <- ggplot(res, aes(x=logFC, y=negLogP)) +
  geom_point(aes(color=negLogP, shape=highlight), alpha=0.85, size=2) +
  scale_color_viridis(option = "D") +
  scale_shape_manual(values = c(19,17)) +
  theme_minimal(base_size = 14) +
  xlab('Log2 fold change (trend)') + ylab('-log10(p-value)') +
  ggtitle('Limma trend: Volcano plot') +
  theme(legend.position='right')

# label top hits by adj.P.Val
top_hits <- head(res[order(res$adj.P.Val), ], n = TOP_LABEL)
# increase repel parameters to avoid overlaps
p <- p + geom_text_repel(data=top_hits, aes(label=protein), size=3.6, max.overlaps = 200, box.padding = 0.3)

# add a horizontal guide for p = 0.001
p <- p + geom_hline(yintercept = -log10(0.001), linetype = 'dashed', color = 'grey50')

# save
png(volcano_png, width=9, height=6, units='in', res=300)
print(p)
dev.off()
cat('Wrote volcano PNG:', volcano_png, '\n')

# Barplot: top absolute logFC (with SE from t-stat)
res$absFC <- abs(res$logFC)
bar_df <- head(res[order(-res$absFC), ], n = TOPN_BAR)
# estimate SE from t-stat (SE = logFC / t) when t != 0
bar_df$SE <- with(bar_df, ifelse(t != 0, logFC / t, NA))
bar_df$protein <- factor(bar_df$protein, levels = rev(bar_df$protein))

bar_png <- file.path(OUT_DIR, 'limma_top_effects_bar_pub.png')
bar_pdf <- file.path(OUT_DIR, 'limma_top_effects_bar_pub.pdf')

bp <- ggplot(bar_df, aes(x=protein, y=logFC, fill=logFC>0)) +
  geom_bar(stat='identity', width=0.7) +
  geom_errorbar(aes(ymin=logFC - SE, ymax=logFC + SE), width=0.3, na.rm=TRUE) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  xlab('Protein') + ylab('Log2 fold change (trend)') +
  ggtitle(sprintf('Top %d proteins by |logFC|', TOPN_BAR)) +
  scale_fill_manual(values = viridis(2, option = "D")) +
  theme(legend.position='none', axis.text.y = element_text(size=9))

png(bar_png, width=8, height=6, units='in', res=300)
print(bp)
dev.off()

pdf(bar_pdf, width=8, height=6)
print(bp)
dev.off()
cat('Wrote barplots:', bar_png, bar_pdf, '\n')
