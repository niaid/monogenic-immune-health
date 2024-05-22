# Utility functions for plotting the Figure 2 DE Heatmaps
# Note that this library has notable changes frequently that may cause problems version to version. Version 2.3.3 is used on my local.
library(ComplexHeatmap)
# This library is important for making the custom dendrogram for our heatmaps
library(dendextend)
# This library is important for the DE colors
library(circlize)
# This utility functions hold the AI/PID/TERT-TERC condition super-types we have been using
source('scripts/util/Groups/groups.R')

util.plot.de.heatmap = function(stats.healthy,
                                stats.all,
                                name, condition.count,
                                pnt.sz = 20, fnt.sz = 12,
                                stat = 't', min.count = 5) {
  ## This utility function plots side-by-side heatmaps for the versus-healthy and versus-all results for each condition,
  ## as well as small heatmaps underneath each of the main ones that display the super-type signatures
  ## Inputs:
  ## stats.healthy - a list of the statistic matrices associated with the versus-healthy DE testing
  ## stats.all - a list of the statistic matrices associated with the versus-all DE testing
  ## name - the name we wish to give the plot
  ## condition.count - a named character vector that maps from conditions to the number of samples associated with that conditon
  ## in the data type
  ## pnt.sz - the size to use for the *'s marking significance in the heatmap
  ## fnt.sz - the font size to use for the feature name labels on the heatmap
  ## stat - the name of the stastic matrix to plot (e.g. 't', 'effect.size', etc.). Should be one of the names of the elements of
  ## stats.healthy and stats.all
  ## min.count - The minimum number of subjects that should be associated with a condition in the data type for that condition to be plotted
  
  ## Get the vectors of conditions in each condition supertype
  AI = util.get_ai()
  PID = util.get_pid()
  TERT.TERC = util.get_tert_terc()
  
  ## Add the AI, PID, and TERT.TERC super-type condition counts
  condition.count['AI'] = sum(condition.count[names(condition.count) %in% AI])
  condition.count['PID'] = sum(condition.count[names(condition.count) %in% PID])
  condition.count['TERT.TERC'] = sum(condition.count[names(condition.count) %in% TERT.TERC])
  
  ## Create a named vector that associates each condition to a string '(n)', where n is the number of subjects associated with that conditon 
  condition_count_str = as.character(condition.count)
  condition_count_str = paste0(' (', condition_count_str, ')')
  names(condition_count_str) = names(condition.count)
  
  ## Get the statistic matrix (we call this effect.sizes but it could in fact be the t-statistic, cohen's d, etc.)
  effect.sizes = t(stats.healthy[[stat]])
  
  ## Decide which conditions pass the minimun count threshold
  conditions.to.keep = condition.count[rownames(effect.sizes)] >= min.count
  ## Append the condition counts to their corresponding conditions
  rownames(effect.sizes) = paste0(rownames(effect.sizes), condition_count_str[rownames(effect.sizes)])
  ## Get the significance matrix associated with the DE testing
  ps1 = t(stats.healthy$adj.P.Val)
  ## Append the condition counts to this matrix as well
  rownames(ps1) = paste0(rownames(ps1, condition_count_str[rownames(ps1)]))
  
  ## Subset the statistic matrix and the pvalue matrix to the desired conditions
  effect.sizes = effect.sizes[conditions.to.keep, , drop = FALSE]
  ps1 = ps1[conditions.to.keep, , drop = FALSE]
  
  ## Find if each row is associated with the condition super-types
  group_rows = rownames(effect.sizes) %in% c(paste0('AI',condition_count_str['AI']), 
                                             paste0('PID', condition_count_str['PID']),
                                             paste0('TERT.TERC', condition_count_str['TERT.TERC']))
  
  ## If there is more than one feature in the matrix
  if(ncol(effect.sizes) > 1) {
    ## Create a dendrogram for this matrix, not including the group rows
    mat = effect.sizes[! group_rows, , drop = FALSE]
    hc = as.dendrogram(hclust(dist(t(mat), method="euclidean"), method = "complete"))
  } else {
    hc = FALSE
  }
  
  ## Create a row grouping to separate the condition rows from condition-supertype rows
  groups = factor(group_rows)
  ## Change the TRUE/FALSE levels of groups to an empty string to avoid labeling the heatmaps separately
  levels(groups) = c(' ', '')
  
  ## Only show column names if there are fewer than 60 features
  show_colnames = ncol(effect.sizes) < 60
  
  ## Create the first heatmap using the versus-healthy option
  H1 = Heatmap(effect.sizes, name = 'healthy', cluster_columns = hc, cluster_rows = T, 
               show_column_names = show_colnames, show_row_names = F,
               split = groups,
               column_title = paste0(name, ' (versus healthy)'),
               column_names_gp = gpar(fontsize = fnt.sz),
               show_heatmap_legend = F,
               col = colorRamp2(c(-4,0,4), c("blue", "white", "red")),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(ps1[i, j] < 0.001) {
                   grid.text("***", x, y, gp = gpar(fontsize = pnt.sz)) ### Use three *'s for p-values less than .001
                 } else if(ps1[i, j] < 0.01) {
                   grid.text("**", x, y, gp = gpar(fontsize = pnt.sz)) ### Use two *'s for p-values less than .01
                 } else if(ps1[i, j] < 0.05) {
                   grid.text("*", x, y, gp = gpar(fontsize = pnt.sz)) ### Use one * for p-values less than .05
                 }
               })
  
  ## Do the same for the versus-healthy matrix
  effect.sizes = t(stats.all[[stat]])
  conditions.to.keep = condition.count[rownames(effect.sizes)] >= min.count
  rownames(effect.sizes) = paste0(rownames(effect.sizes), condition_count_str[rownames(effect.sizes)])
  ## Note that due to lazy evaluation of functions in R, this pvalue matrix must be named diffently from ps1, or the previous heatmap's cell_fun will use the newer matrix
  ps2 = t(stats.all$adj.P.Val)
  rownames(ps2) = paste0(rownames(ps2, condition_count_str[rownames(ps2)]))
  
  effect.sizes = effect.sizes[conditions.to.keep, , drop = FALSE]
  ps2 = ps2[conditions.to.keep, , drop = FALSE]
  
  H2 = Heatmap(effect.sizes, name = 'H2', column_order = column_order(H1), 
               cluster_columns = F, cluster_rows = F,
               show_column_names = show_colnames,
               column_title = paste0(name, ' (versus all conditions)'),
               column_names_gp = gpar(fontsize = fnt.sz),
               row_names_max_width = max_text_width(rownames(effect.sizes)),
               col = colorRamp2(c(-4,0,4), c("blue", "white", "red")),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(ps2[i, j] < 0.001) {
                   grid.text("***", x, y, gp = gpar(fontsize = pnt.sz))
                 } else if(ps2[i, j] < 0.01) {
                   grid.text("**", x, y, gp = gpar(fontsize = pnt.sz))
                 } else if(ps2[i, j] < 0.05) {
                   grid.text("*", x, y, gp = gpar(fontsize = pnt.sz))
                 }
               },
               heatmap_legend_param = list(color_bar = "continuous",
                                           legend_width = unit(3,"cm"),
                                           legend_direction = "vertical"))
  
  # Combine the heatmaps
  htlist = H1 + H2
  # And draw then
  draw(htlist)
  return('Complete')
}

util.plot.de.heatmap.2 = function(stats.healthy,
                                  stats.all,
                                  name, pnt.sz = 20, fnt.sz = 12) {
  ## This function is similar to the last, except that it does not separate the condition super-type rows from the
  ## condition rows. 
  
  effect.sizes = t(stats.healthy$t)
  ps1 = t(stats.healthy$adj.P.Val)
  
  if(ncol(effect.sizes) > 1) {
    mat = effect.sizes[! rownames(effect.sizes) %in% c('AI', 'PID', 'TERT.TERC'), ]
    hc = as.dendrogram(hclust(dist(t(mat), method="euclidean"), method = "complete"))
  } else {
    hc = FALSE
  }
  
  show_colnames = ncol(effect.sizes) < 60
  
  H1 = Heatmap(effect.sizes, name = 'healthy', cluster_columns = hc, cluster_rows = T, 
               show_column_names = show_colnames, show_row_names = F,
               column_title = paste0(name, ' (versus healthy)'),
               column_names_gp = gpar(fontsize = fnt.sz),
               show_heatmap_legend = F,
               col = colorRamp2(c(-4,0,4), c("blue", "white", "red")),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(ps1[i, j] < 0.001) {
                   grid.text("***", x, y, gp = gpar(fontsize = pnt.sz))
                 } else if(ps1[i, j] < 0.01) {
                   grid.text("**", x, y, gp = gpar(fontsize = pnt.sz))
                 } else if(ps1[i, j] < 0.05) {
                   grid.text("*", x, y, gp = gpar(fontsize = pnt.sz))
                 }
               })
  
  effect.sizes = t(stats.all$t)
  ps2 = t(stats.all$adj.P.Val)
  H2 = Heatmap(effect.sizes, name = '', column_order = column_order(H1), 
               cluster_columns = F, cluster_rows = F,
               show_column_names = show_colnames,
               column_title = paste0(name, ' (versus all conditions)'),
               column_names_gp = gpar(fontsize = fnt.sz),
               row_names_max_width = max_text_width(rownames(effect.sizes)),
               col = colorRamp2(c(-4,0,4), c("blue", "white", "red")),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if(ps2[i, j] < 0.001) {
                   grid.text("***", x, y, gp = gpar(fontsize = pnt.sz))
                 } else if(ps2[i, j] < 0.01) {
                   grid.text("**", x, y, gp = gpar(fontsize = pnt.sz))
                 } else if(ps2[i, j] < 0.05) {
                   grid.text("*", x, y, gp = gpar(fontsize = pnt.sz))
                 }
               },
               heatmap_legend_param = list(color_bar = "continuous",
                                           legend_width = unit(3,"cm"),
                                           legend_direction = "vertical"))
  
  htlist = H1 + H2
  draw(htlist)
  return('Complete')
}

