library(tidyverse)
library(variancePartition)
library(knitr)

setwd("/hpcdata/sg/sg_data/PROJECTS/Monogenic_Project")
source("scripts_nick/util/Enrichment/plot_camera_enrich.R")


#ifn.grepnames <- c("M150", "M75", "INTERFERON_ALPHA_BETA_SIGNALING", "INTERFERON_GAMMA_SIGNALING", "INTERFERON_SIGNALING", "IFN.I.Dcact")
#ifn.genesets <-sapply(ifn.grepnames, grep, enrich.dat.list$XCGD$geneset, value = TRUE, ignore.case = T)
ENRICHMENT.DAT.IN.PATH <- "Data/Microarray/analysis_output/differential_expression/enrichment/cameraPR_enrichment_list.RDS"

ALL.FIG.OUT.PATH <- "Misc/plots/figures_for_rachel/interferon_enrichments/array_de_enrich.pdf"
#IFN.FIG.OUT.PATH <- "Misc/plots/figures_for_rachel/interferon_enrichments/array_de_enrich_ifn_only.pdf"
dir.create(dirname(FIG.OUT.PATH), recursive = T)

enrich.dat.list <- readRDS(ENRICHMENT.DAT.IN.PATH)

pdf(ALL.FIG.OUT.PATH)
for(nm in names(enrich.dat.list)){
        plot_enrich_both_directions(enrich.dat.list[[nm]], main = nm)
}
dev.off()

#pdf(IFN.FIG.OUT.PATH)
#for(nm in names(enrich.dat.list)){
#        dat <- enrich.dat.list[[nm]]
#        dat <- dat[dat$geneset %in% ifn.genesets, ]
#        plot_enrich_both_single_plot(dat, main = nm)
#}
#dev.off()
