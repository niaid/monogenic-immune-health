library(tidyverse)

setwd("/hpcdata/sg/sg_data/PROJECTS/Monogenic_Project/")

soma <- readRDS("Data/Somalogic/data_analysis_ready/analysis_ready_sample_level_training_somalogic.rds")
somamers_filtered <- data.table::fread("Data/Somalogic/raw/v1/Cal_QC_CHI_Somamers.txt")
nrow(soma)

nrow(soma2)


somamers_orig <- readr::read_tsv("Somalogic/Tsang/CHI-16-017/DATA/CHI-16-017-SET003.20161222.adat_Somamers.txt", col_names =FALSE) %>% as.data.frame()  %>% as.matrix()

somamers_orig <- somamers_orig %>%
        t()
somamers_orig <- somamers_orig %>%
        #`rownames<-`(.[,1]) %>%
        `colnames<-`(.[1,]) %>%
        .[-1, ]

somamers_orig <- somamers_orig %>% as.data.frame(stringsAsFactors = FALSE)

table(somamers_orig$Type)
table(somamers_filtered$Type)
somamers_orig %>% filter(SomaId %in% somamers_filtered$SomaId) %>%
        pull(Type) %>%
        table()

grep("HPV", somamers_orig$Target, value = T, ignore.case = T)
grep("EGFRvIII", somamers_orig$Target, value = T, ignore.case = T)
grep("P05186", somamers_orig$UniProt, value = T, ignore.case = T)
somamers_orig$UniProt[somamers_orig$Target == "EGFRvIII"]

length(unique(somamers_filtered$Target))
length(unique(somamers_filtered$UniProt))



