library(tidyverse)

load("/Volumes/SG_DATA/PROJECTS/Monogenic_Project/Metadata/monogenic.de-identified.metadata.RData")

dat <- monogenic.all.assays %>%
  filter(analysis_group == "Discovery")

# visit_count <- dat %>%
#   group_by(patient_id) %>%
#   summarise()
# 
# rep_pats <- 

summ <- dat %>%
  #filter(patient_id %in% rep_pats) %>%
  group_by(patient_id) %>%
  summarise(min_day = min(visit_date), max_day = max(visit_date),
            n_visit = length(unique(visit_id))) %>%
  mutate(delta_days = max_day - min_day) %>%
  filter(n_visit > 1) %>%
  ungroup()

hist(as.numeric(summ$delta_days), breaks = 100)

