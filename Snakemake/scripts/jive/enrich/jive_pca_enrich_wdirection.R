# Computes wilcoxon based geneset enrichments for JIVE PC scores 
# for joint and individual PC's.
# Iterates through various geneset databases.
# All output is saved in a sinlge large dataframe that is used for figure generation.


suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(parallel)
  library(Biobase)
  library(BiocGenerics)
})

#Source Scripts ----------------------------------------------------------
#source("scripts/util/Processing/averageRepeatSamples.R")
source("scripts/util/Enrichment/camera.R")
source("scripts/util/JIVE/prcomp_list_varfilter.R")

#Set paths ---------------------------------------------------------------
#Jive
JIVE.PATH <- snakemake@input[["jive"]]#"Integration_output/jive/subject/jive.rds"
JIVE.PC.PATH <- snakemake@input[["jive_pcs"]]#"Integration_output/jive/subject/prcomp_list.rds"

#Data
ARRAY.ESET.PATH <- snakemake@input[["array_eset"]]#"Data/Microarray/data_analysis_ready/eset_batch_training_subject.rds"
SOMA.ESET.PATH <- snakemake@input[["soma_eset"]]#"Data/Somalogic/data_analysis_ready/analysis_ready_sample_level_training_somalogic.rds"

#Geneset paths
TISSUE.IN.PATH <- snakemake@input[["tissue_genesets"]]#"Gene_sets/processed/tissue_gene_sets.RDS"
COMBINED.GENESETS.IN.PATH <- snakemake@input[["genesets"]]#"Gene_sets/processed/combined_gene_sets.RDS"

#Out path
ENRICHMENT.DAT.OUT.PATH <- snakemake@output[[1]]#"Integration_output/jive/subject/pc_enrich_dat_camera.rds"

#Load Jive ---------------------------------------------------------------
#need for the correlation with the array and soma data
jive <- readRDS(JIVE.PATH)

#Load Jive pcs------------------------------------------------------------
prcomp.list <- readRDS(JIVE.PC.PATH)
pdat <- prcomp.list$pdat

#remove unneccessary PC's ----------------------------------------------
prcomp.list <- prcomp_list_varfilter(prcomp.list, .03)
pca.list <- map(prcomp.list, "x") # get the PC's

#load data that has not been filtered for stable features ----------------
array.eset <- readRDS(ARRAY.ESET.PATH)
soma.eset <- readRDS(SOMA.ESET.PATH)

#Load Gene sets ----------------------------------------------------------
genesetLL <- readRDS(COMBINED.GENESETS.IN.PATH)
tissue <- readRDS(TISSUE.IN.PATH)

#create geneset list of lists --------------------------------------------
genesetLL <- c(genesetLL, list(tissue = tissue$general))

#Change the somalogic names so that they are genes not proteins ----------
soma.fdata <- featureData(soma.eset)@data
soma.data <- t(jive$data$soma)
array.data <- t(jive$data$array)
colnames(soma.data) <- 
  soma.fdata$EntrezGeneSymbol[match(colnames(soma.data), soma.fdata$Target)]
in.data.list <- list(array = array.data, soma = soma.data)

#Run enrichment -----------------------------------------------------------
enrich.dat <- lapply(pca.list, function(pca){
  lapply(genesetLL, function(geneset.list){
    lapply(in.data.list, function(in.data){
      cormat <- get_cormat(pca, in.data)
      
      cameraPR_cor(cormat, geneset.list, 
                   min.geneset.size = 3, use.ranks = FALSE, 
                   abs.cor = F) %>% 
        bind_rows(.id = "PC")
    }) %>% bind_rows(.id = "in.data")
  }) %>% bind_rows(.id = "geneset.db")
}) %>% bind_rows(.id = "pca.data")

#Save ----------------------------------------------------------------------
saveRDS(enrich.dat, ENRICHMENT.DAT.OUT.PATH)

