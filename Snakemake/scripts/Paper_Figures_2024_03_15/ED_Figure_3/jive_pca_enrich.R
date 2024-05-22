library(tidyverse)

DAT.IN.PATH <- snakemake@input[[1]]#"Integration_output/jive/subject/pc_enrich_dat_camera.rds"

FIG.OUT.PATH <- snakemake@output[[1]]#"Paper_1_Figures/Supplemental_Figure_3/jive_pc_enrichments_joint_down_pc1.pdf"

all.dat <- readRDS(DAT.IN.PATH)
all.dat <- all.dat %>% filter(geneset.db != "tiss.general")
all.dat$col <- as.numeric(as.factor(all.dat$geneset.db))

plot_single <- function(dat, 
                        n.genesets = 20,
                        main = ""
                        ){
  #filter to the top ranked by FDR for that comparison
  selection2 <- order(dat[["FDR_all_db"]])[seq_len(n.genesets)]
  dat <- dat[rev(selection2), ]

  dat <- dat %>%
          mutate(geneset = paste(geneset.db, geneset, sep = "_"))

  #geneset_vec <- strsplit(dat$geneset, split = "_")

  #nwords_in_geneset_name <- sapply(geneset_vec, length)

  #for(i in seq_along(geneset_vec)){
  #  #print(nwords_in_geneset_name[[i]])
  #  if(nwords_in_geneset_name[[i]] > 6){
  #        geneset_vec[[i]] <- append(geneset_vec[[i]], "\n", after = 6)
  #      }
  #  geneset_vec[[i]] <- paste(geneset_vec[[i]], collapse = "_")
  #}

  #dat$geneset <- geneset_vec

  #if(nwords_in_geneset_name)

  geneset_levels <- dat$geneset[order(dat[["FDR_all_db"]])]

  dat <- dat %>%
          mutate(geneset = factor(geneset, levels = geneset_levels))

  p <- ggplot(dat, aes(y = -log10(FDR_all_db), x = geneset)) + 
          geom_col() +
          theme_bw() +
          xlab("") +
          coord_flip() +
          ggtitle(main)
  return(p)
}


all.dat <- all.dat %>%
        group_by(pca.data, in.data, PC, Direction) %>%
        mutate(FDR_all_db = p.adjust(PValue, method = "fdr"))


plot_list <- list()
i <- 0
for(indata in c("array", "soma")){
  for(PC. in c("PC1")){
    for(direct in c("Up", "Down")){
     
      p <- all.dat %>%
              filter(pca.data == "joint", in.data == indata, 
                     PC == PC., Direction == direct) %>% 
              plot_single(main = paste("joint", indata, PC., direct), n.genesets  = 20)
    
      i <- i + 1
      plot_list[[i]] <- p
    
    }
  }
}
library(cowplot)
p_all <- plot_grid(plotlist = plot_list, ncol =2, nrow = 1)

pdf(FIG.OUT.PATH, 
    height = 3.5, width = 11)
print(p_all)

dev.off()

