library(GEOquery)
library(Biobase)

# Set Globals
CGPS.IN.PATH = 'Reference/jamboree/processed/jamboree_cgps.RDS'
GEO.OUT.DIR = 'Reference/jamboree/raw/series_matrices'

cgps <- readRDS(CGPS.IN.PATH)

#save one object that has the series matrix saved with the cgp -------------------------
for(nm1 in names(cgps)){
  for(nm2 in names(cgps[[nm1]])){
    geo.id <- cgps[[nm1]][[nm2]]$study.info$study
    geo.platform <- cgps[[nm1]][[nm2]]$study.info$platform
    geo.query.object <- getGEO(geo.id, destdir = GEO.OUT.DIR)
    #some geoquery objects have multiple platforms. I pull the one that is listed in our CGP object
    if(length(geo.query.object) == 1){
      #if there are not multiple platforms, just select the single one
      keep.platform.index <- 1
    }else{
      keep.platform.index <- which(grepl(geo.platform, names(geo.query.object)))
    }
    series.matrix <- phenoData(geo.query.object[[keep.platform.index]])@data
    out.name <- paste(geo.id, geo.platform, 'series_matrix.txt', sep = '_')
    out.path <- file.path(GEO.OUT.DIR, out.name)
    write.table(series.matrix, out.path, row.names = TRUE, col.names = TRUE, quote = FALSE, comment.char = '', sep = '\t')
  }
}
