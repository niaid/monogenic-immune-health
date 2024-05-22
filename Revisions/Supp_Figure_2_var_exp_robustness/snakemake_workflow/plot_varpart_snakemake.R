library(variancePartition)
library(SummarizedExperiment)

# Input Paths -----------------------------------------------------
AGE.FIRST.PAX <- snakemake@params[["afp"]]
print(AGE.FIRST.PAX)
VARPART.IN.PATH <- snakemake@input[["varpart"]]
VST.IN.PATH <- snakemake@input[["vst"]]

# Output Paths
FIG.OUT.PATH <- snakemake@output[[1]]

# Read in data ----------------------------------------------------
vp <- readRDS(VARPART.IN.PATH)
vst <- readRDS(VST.IN.PATH)

pdf(FIG.OUT.PATH)
#plot with all genes
plotVarPart(vp, main = "All genes")

#Filtering to the top 20% of variable genes

# subset to particular age first pax
vst <- vst[, as.character(vst$age_first_pax) %in% AGE.FIRST.PAX]

vst.mat <- assays(vst)$vst
sds <- apply(vst.mat, 1, sd)
var.genes <- names(sds)[sds > quantile(sds, .8, na.rm = TRUE)]
vp.sub <- subset(vp, rownames(vp) %in% var.genes)

#plot with variable genes
plotVarPart(vp.sub, main = "top 20% most variable genes")

dev.off()


