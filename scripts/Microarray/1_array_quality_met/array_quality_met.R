library(arrayQualityMetrics)

#config---------------------------------------------
#input
ESET.PATH <- "Data/Microarray/raw/eset_rma_with_pData.rds"
#output
REPORT.DIR <- "Data/Microarray/Array_quality_metrics"
OUT.AQM.PATH <- file.path(REPORT.DIR, "aqm_table.rds")

#load eset-------------------------------------------
eset <- readRDS(ESET.PATH)

#run arrayQualityMetrics-----------------------------
aqm <- arrayQualityMetrics(expressionset = eset,
                    outdir = REPORT.DIR,
                    force = TRUE, spatial = FALSE)

aqm_table <- aqm$arrayTable

colnames(aqm_table)[colnames(aqm_table) == '<a href=\"#hm\">*1</a>'] <- "abs_mean_distance"
colnames(aqm_table)[colnames(aqm_table) == '<a href=\"#box\">*2</a>'] <- "KS_Statistic_on_signalintensity"
colnames(aqm_table)[colnames(aqm_table) == '<a href=\"#ma\">*3</a>'] <- "D_statistic_on_MAplot"

saveRDS(aqm_table, OUT.AQM.PATH)







