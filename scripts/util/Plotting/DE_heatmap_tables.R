library(ComplexHeatmap)
library(dendextend)

util.plot.de.heatmap = function(stats.healthy,
                                stats.all,
                                name, pnt.sz = 20, fnt.sz = 12, stats = 'cohens.d') {
  
  
  effect.sizes = t(stats.healthy[[stats]])
  effect.sizes = effect.sizes[! rownames(effect.sizes) %in% c('AI', 'PID', 'TERT.TERC'), ]
  ps1 = t(stats.healthy$adj.P.Val)
  ps1 = ps1[! rownames(ps1) %in% c('AI', 'PID', 'TERT.TERC'), ]

  hc.col = hclust(dist(t(effect.sizes), method="euclidean"), method = "complete")
  hc.row = hclust(dist(effect.sizes, method="euclidean"), method = "complete")
  
  show_colnames = ncol(effect.sizes) < 60
  
  H1 = Heatmap(
    effect.sizes, name = 'healthy', cluster_columns = hc.col, cluster_rows = hc.row,
    show_column_names = show_colnames, show_row_names = T,
    column_title = paste0(name, ' (versus healthy)'),
    column_names_gp = gpar(fontsize = fnt.sz),
    show_heatmap_legend = T,
    col = colorRamp2(c(-4,0,4), c("blue", "white", "red"))
  )
  
  table.rows = rownames(effect.sizes)[hc.row$order]
  table.cols = colnames(effect.sizes)[hc.col$order]
  effect.sizes.table = effect.sizes[table.rows, table.cols]
  
  draw(H1)
  return(effect.sizes.table)
}
