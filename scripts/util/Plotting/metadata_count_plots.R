suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

metadata_plot_multi <- function(expressionset, traits.of.interest, plot.type = c("bar", "dot")){
  #returns a list of ggplot objects than will be plotted
  bar.plots <- lapply(traits.of.interest, function(trait){
    if(plot.type == "bar"){
      single.plot <- metadata_bar_plot(expressionset, trait)
    }else if(plot.type == "dot"){
      single.plot <- metadata_dot_plot(expressionset, trait)
    }else{
      print("Invalid plot type. Must be 'bar' or 'dot'")
      return(NULL)
    }
    return(single.plot)
  }
  )
  return(bar.plots)
}

metadata_bar_plot <- function(expressionset, trait){
  
  ## traits of interest is character vector
  pdat <- pData(eset) %>% select(trait)
  
  g <- ggplot(pdat, aes_string(x = trait)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(trait)
  
  # print(trait)
  # print(class(pdat[[trait]]))
  
  if(class(pdat[[trait]])== "factor"){
    gg <- g + geom_bar()
  }else if(class(pdat[[trait]])== "numeric"){
    gg <- g + geom_histogram()
  }
  
  return(gg)
}

metadata_dot_plot <- function(expressionset, trait){
  ## Plots counts of combinations of patient id, condition, and single traits of interest
  
  pdat <- pData(expressionset) %>% select(c(trait, "patient_id", "condition"))
  # print(str(pdat))
  
  g <- ggplot(pdat, aes_string(x = trait)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(trait)
  
  # print(class(pdat[[trait]]))
  # print(trait)
  
  ggplot(pdat, aes(x = patient_id, y = condition)) + 
    geom_count(aes_string(color = trait)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
}
