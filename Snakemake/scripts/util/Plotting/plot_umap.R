library(umap)
library(ggplot2)

plot_umap <- function(mat, color = NULL, label = NULL){
  umap.out <- umap(mat, method = "umap-learn")$layout
  qplot(umap.out[,1], umap.out[,2], color = color, geom = "text", label = label)
}
