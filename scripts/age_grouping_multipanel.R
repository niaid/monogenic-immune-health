# scripts/age_grouping_multipanel.R
# Combine heatmap + volcano + barplot PNGs into a single multipanel PDF

required_pkgs <- c('png','grid','gridExtra')
for(pk in required_pkgs){
  if(!requireNamespace(pk, quietly = TRUE)){
    install.packages(pk, repos = 'https://cran.rstudio.com/')
  }
}

library(png)
library(grid)
library(gridExtra)

OUT_DIR <- 'outputs'
heat_png <- file.path(OUT_DIR, 'age_diff_heatmap_limma_pub.png')
volcano_png <- file.path(OUT_DIR, 'limma_volcano_pub.png')
bar_png <- file.path(OUT_DIR, 'limma_top_effects_bar_pub.png')

missing <- c(heat_png, volcano_png, bar_png)[!file.exists(c(heat_png, volcano_png, bar_png))]
if(length(missing)>0){
  stop(sprintf('Missing expected images: %s', paste(missing, collapse=', ')))
}

# read images
heat_img <- readPNG(heat_png)
volcano_img <- readPNG(volcano_png)
bar_img <- readPNG(bar_png)

# convert to grobs
g_heat <- rasterGrob(heat_img, interpolate=TRUE)
g_volcano <- rasterGrob(volcano_img, interpolate=TRUE)
g_bar <- rasterGrob(bar_img, interpolate=TRUE)

# arrange: heatmap left (two-thirds height), volcano + bar stacked right
# For simplicity, use a 2-column layout: heatmap (col1 full height), right column stacked volcano then bar

out_pdf <- file.path(OUT_DIR, 'age_analysis_multipanel_pub.pdf')
# choose page size (wide): width 14in x height 8.5in
pdf(out_pdf, width = 14, height = 8.5)

# layout matrix: 2 cols, 2 rows, heatmap spans rows 1:2 col1
layout_mat <- rbind(c(1,2), c(1,3))

grid.arrange(grobs = list(g_heat, g_volcano, g_bar), layout_matrix = layout_mat,
             widths = unit(c(2/3,1/3), 'npc'))

dev.off()
cat('Wrote multipanel PDF:', out_pdf, '\n')
