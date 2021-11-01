
library(tidyverse)

load("../../Inputs/Metadata/monogenic.de-identified.metadata.RData")

conditions = monogenic.all.assays[monogenic.all.assays$analysis_group %in% c('Discovery'), 'condition'] %>% unique

## Remove any unknown samples from these conditions
unknowns = c('Unknown', 'U_PI3K', 'U_Kastner', 'U_STAT1', 'U_CTLA4', 'U_Telomere')
conditions = setdiff(conditions, unknowns)

summ <- monogenic.all.assays %>%
        filter(assay_type %in% c("Microarray", "Somalogic")| patient_id %in% c(42, 59, 78 171)) %>%
        group_by(patient_id, visit_id) %>%
        summarise(visit_dates = paste(unique(visit_date), collapse = ";"),
        blood_draw_dates = paste(unique(blood_draw_date), collapse = ";"),
        diff_blood_draw = max(as.numeric(blood_draw_date)) - min(as.numeric(blood_draw_date))) %>%
        ungroup()

summ$diff_blood_draw


table(summ$diff_blood_draw > 0)
table(summ$diff_blood_draw > .5)

