library(tidyverse)
library(cowplot)
library(ggpubr)

ALL.SUBJECT.PATH <- snakemake@input[["all_subj"]]#"Integration_output/jive/subject/prcomp_list.rds"
HEALTHY.ONLY.PATH <- snakemake@input[["healthy_only"]]#"Integration_output/jive/subject_onlyHealthy/prcomp_list.rds"
NO.HEALTHY.PATH <- snakemake@input[["no_healthy"]]#"Integration_output/jive/subject_noHealthy/prcomp_list.rds"

PLOT_OUT_PATH <- snakemake@output[[1]]#"Paper_1_Figures/Figure_3/jive_pc_comparison.pdf"

path.list <- list(`All Subjects` = ALL.SUBJECT.PATH, 
                  `Only Healthy` = HEALTHY.ONLY.PATH, 
                  `No Healthy` = NO.HEALTHY.PATH)

prcompLL <- lapply(path.list, readRDS)

joint.list <- lapply(prcompLL, function(prcompL){
  pca <- prcompL[["joint"]]
  
  pca[["x"]]
})

names(joint.list) <- gsub("No Healthy", "Excluding Healthy", names(joint.list))
names(joint.list) <- gsub("All Subjects", "Including All Subjects", names(joint.list))

#This function takes of the jive categories and plots the desired pc
# Scatterplot with correlation p values

single_plot <- function(jive.cat1, jive.cat2, PC){
  x <- joint.list[[jive.cat1]][, PC]
  y <- joint.list[[jive.cat2]][, PC]
  
  intersecting.pats <- intersect(names(x), names(y))
  x <- x[intersecting.pats]
  y <- y[intersecting.pats]
  
  #flip direction if anticorrelated- directions of PC's are arbitrary
  if(cor(x, y) <= 0){
    y <- -y
  }
  
  dat <-data.frame(x =x, y = y)
  p <- ggplot(dat, aes(x =x, y=y))+
    geom_point()+
    xlab(jive.cat1) +
    ylab(jive.cat2) + 
    stat_smooth(method = "lm", se = FALSE) +
    ggtitle(PC) + 
    theme_bw() +
    stat_cor()
  return(p)
}

#Creat PC1-3 plot objects to be printed later

pc1.plotlist <- list(single_plot("Including All Subjects", "Excluding Healthy", "PC1"),
                     single_plot("Including All Subjects", "Only Healthy", "PC1"))
  



#final plot used in figure

p <- plot_grid(plotlist = pc1.plotlist, ncol = 2)
ggsave(plot = p, PLOT_OUT_PATH, height = 4)
