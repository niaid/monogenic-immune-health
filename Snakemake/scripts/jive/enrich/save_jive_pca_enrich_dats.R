IN.DAT.PATH <- "Integration_output/jive/subject/pc_enrich_dat_camera.rds"
OUT.PATH <- "Integration_output/jive/subject/pc_enrich_dats"

dat <- readRDS(IN.DAT.PATH)
dat <- dat %>% filter(geneset.db != "tiss.general")

if(!dir.exists(OUT.PATH)){
  dir.create(OUT.PATH)
}

save_dat <- function(dat, out.dir.path){
  for(in.data in unique(dat$in.data)){
    path1 <- file.path(out.dir.path, in.data)
    dir.create(path1)
    for(pca.data in unique(dat$pca.data)){
      path2 <- file.path(path1, pca.data)
      dir.create(path2)
      for(PC in unique(dat$PC)){
        path3 <- file.path(path2, PC)
        dir.create(path3)
        for(direction in unique(dat$Direction)){
          final.path <- file.path(path3, direction)
          final.path <- paste0(final.path, ".csv")
          #filter to desired comparison
          selection1 <- 
            dat$in.data == in.data &
            dat$pca.data == pca.data &
            dat$PC == PC &
            dat$Direction == direction
          dat.sub <- dat[selection1, ]
          #Adjust the p values accordingly
          dat.sub$FDR <- p.adjust(dat.sub$PValue, method = "BH")
          #sort by FDR
          dat.sub <- dat.sub[rev(order(dat.sub$FDR)), ]
          #save the data
          print(final.path)
          write.csv(dat.sub, final.path)
        }
      }
    }
  }
}

save_dat(dat, OUT.PATH)
