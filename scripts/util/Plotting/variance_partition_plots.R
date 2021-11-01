suppressPackageStartupMessages({
  library(gplots)
  library(dplyr)
  library(ggplot2)
})

varPartHeatmap <- function(varpart.obj, show.rownames = FALSE){
  if(show.rownames == FALSE){
    labRow <- ""
  }else{
    labRow <- NULL
  }
  
  heatmap.2(as.matrix(data.frame(varpart.obj)), 
            scale = "none", 
            labRow = labRow,
            cexCol = .5, 
            srtCol = 20)
}

module_variance_plot <- function(expressionset, i){
  dat <- pData(expressionset) %>% mutate(activity = exprs(expressionset)[i,])
  #print(colnames(dat))
  gg <- 
    ggplot(dat, aes(x = patient_id, y = activity)) + 
    geom_point(aes(color = patient_id), show.legend = FALSE) + 
    theme(axis.text.x = element_blank()) +
    facet_grid(~condition) + #works on HPC
    #facet_grid(1~condition) + #works on my local
    ggtitle(featureNames(expressionset)[i])
  return(gg)
}

module_variance_plot_multi <- function(expressionset){
  if(nrow(exprs(expressionset)) > 100){
    break("stop- you are trying to plot an expression set with too many features")
  }else{
    var.plots <-lapply(1:nrow(expressionset), function(i){module_variance_plot(expressionset, i)})
  } 
  return(var.plots)
  
}
