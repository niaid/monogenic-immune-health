library(limma)

#Source Scripts ----------------------------------------------------------
source("scripts/util/Enrichment/camera.R")
#Set paths ---------------------------------------------------------------
#Inputs
FIT.IN.PATH <- "Data/Microarray/analysis_output/differential_expression/microarray_features_DE_fit.rds"
COMBINED.GENESETS.IN.PATH <- "Gene_sets/processed/combined_gene_sets.RDS"

#Out path
ENRICHMENT.DAT.OUT.PATH <- "Data/Microarray/analysis_output/differential_expression/enrichment/cameraPR_enrichment_list.rds" 
dir.create(dirname(ENRICHMENT.DAT.OUT.PATH), recursive = T)

#Load data ---------------------------------------------------------------
fit <- readRDS(FIT.IN.PATH)
genesetLL <- readRDS(COMBINED.GENESETS.IN.PATH)

#remove that t-statistics for genes without name
fit <- eBayes(fit)
tstat.dat <- fit$t
tstat.dat <- tstat.dat[!is.na(rownames(tstat.dat)), ]

#concatenate geneset list into single list
names(genesetLL$reactome) <- paste0("reactome_", names(genesetLL$reactome))
names(genesetLL$btms) <- paste0("btm_", names(genesetLL$btms))

geneset.list <- Reduce(c, genesetLL)

IFN.I.Dcact <- c("ATF3", "CCL8", "CXCL10", "DDX58", "DDX60", "DHX58", "EIF2AK2", 
                 "HERC5", "IFI27", "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IRF7",
                 "LAMP3", "MX2", "OAS1", "OAS3", "OASL", "PARP9", "PLSCR1", "PML", 
                 "RSAD2", "SERPING1", "SP100", "TAP1")
geneset.list <- c(geneset.list, list(baseline_IFN.I.Dcact = IFN.I.Dcact))

#Run enrichment -----------------------------------------------------------
universe <- rownames(tstat.dat)
indices <- make_indices(geneset.list, universe, 5)

# keep only the condition ceofficients
tstat.dat <- tstat.dat[, grep("group", colnames(tstat.dat))]
colnames(tstat.dat) <- gsub("group", "", colnames(tstat.dat))

enrich.dat.list <- lapply(colnames(tstat.dat), function(col.name){
  tstat <- tstat.dat[, col.name]
  
  enrich.dat <- cameraPR(tstat, indices, use.ranks = TRUE, sort = FALSE)
  enrich.dat$geneset <- names(indices)
  geneset.db <- sapply(strsplit(enrich.dat$geneset, "_"), `[[`, 1)
  enrich.dat$geneset.db <- geneset.db
  
  return(enrich.dat)
})
names(enrich.dat.list) <- colnames(tstat.dat)

#Save ----------------------------------------------------------------------
saveRDS(enrich.dat.list, ENRICHMENT.DAT.OUT.PATH)

