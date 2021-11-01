#library(tidyverse)

setwd("/hpcdata/sg/sg_data/PROJECTS/Monogenic_Project")

tbnk <- readRDS("Data/TBNK/data_analysis_ready/tbnk_training_sample_level_eset.rds")

fdat <- featureData(tbnk)@data
