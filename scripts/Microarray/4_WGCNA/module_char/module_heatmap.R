rm(list = ls())

library(Biobase)
library(BiocGenerics)
library(data.table)
library(gplots)

#wgcna.dir <- "/Volumes/SG_DATA/users/rachmaninoffn2/monogenic/WGCNA_array/"
wgcna.dir <- "/hpcdata/sg/sg_data/users/rachmaninoffn2/monogenic/WGCNA_array/"

out.fp <- paste0(wgcna.dir, "wgcna_out/module_heatmap")

load(paste0(wgcna.dir,"processed_eset/filtered_wgcna_eset.rdata")) #expressionset
load(paste0(wgcna.dir,"wgcna_out/wgcna_res.rdata")) #wgcna output

#corMat <- cor(t(exprs(eset[order(moduleColors),])))

corMat <- cor(t(exprs(eset[order(moduleColors), ])))

moduleColors <- moduleColors[order(moduleColors)]

png(filename = out.fp, type = "cairo")
heatmap.2(corMat, srtCol = 20, 
          margin=c(5,15),
          #dendrogram= TRUE,
          #Rowv=FALSE, 
          RowSideColors=moduleColors,
          labRow = FALSE, 
          scale = "none"
)

legend("topright",      
       legend = unique(moduleColors),
       col = unique(moduleColors),
       lty= 1,             
       lwd = 5,           
       cex=.7
)

dev.off()