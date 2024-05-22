plot_enrich_single <- function(dat, direction, main){
  par(mar=c(5, 12, 1, 1))
  tmp <- dat %>% 
    filter(Direction == direction) %>% #Select up or down
    arrange(-FDR) %>% #sort by the FDR
    top_n(30, -FDR)
  
  
  barplot(-log10(tmp$FDR), names.arg = tmp$geneset,
          horiz = TRUE, las = 2, cex.names = .35,
          main = main,
          xlab = "-log10(FDR)")
}

plot_enrich_both_single_plot <- function(dat,  main){
  par(mar=c(5, 12, 1, 1))
  tmp <- dat %>% 
    arrange(-FDR) %>% #sort by the FDR
    mutate(FDR = FDR * ifelse(Direction == "Up", 1, -1))
  
  barplot(-log10(tmp$FDR), names.arg = tmp$geneset,
          horiz = TRUE, las = 2, cex.names = .35,
          main = main,
          xlab = "-log10(FDR)")
}


plot_enrich_both_directions <- function(dat, main){
  par(mar=c(3, 12, 1, 1))
  
  
  for(direction in c("Up", "Down")){
    tmp <- dat %>% 
      filter(Direction == direction) %>% #Select up or down
      arrange(-FDR) %>% #sort by the FDR
      top_n(30, -FDR)
    
    
    barplot(-log10(tmp$FDR), names.arg = tmp$geneset,
            horiz = TRUE, las = 2, cex.names = .35,
            main = paste(main, direction))
  }
}
