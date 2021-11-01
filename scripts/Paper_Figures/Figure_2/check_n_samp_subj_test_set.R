library(tidyverse)
library(Biobase)


test <- readRDS("/hpcdata/sg/sg_data/PROJECTS/Monogenic_Project/Data/Somalogic/data_analysis_ready/analysis_ready_sample_level_testing_somalogic.rds")

meta <- pData(test)

summ <- meta %>%
        group_by(condition) %>%
        summarise(n_subj = length(unique(patient_id)),
                  n_samp = n())

as.data.frame(summ) %>%
        arrange(-n_samp)
