suppressPackageStartupMessages({
  library(Biobase)
  library(BiocGenerics)
  library(dplyr)
})


subsetFirstRoutineVisit <- function(expressionset){
  #select routine(not sick) visits
  expressionset <- expressionset[, expressionset$visit_type == "Routine"]
  
  #find first visit
  first.visits <- sapply(unique(expressionset$patient_id), function(patient){
    dat <- pData(expressionset) %>% filter(patient_id == patient)
    first.visit.date <- min(as.numeric(dat$visit_date))
    first.visit.id <- dat %>% 
      filter(visit_date == first.visit.date) %>%
      select(visit_id)
    return(first.visit.id)
  })
  
  #subset to just first visits for each patient
  expressionset <- expressionset[, expressionset$visit_id %in% first.visits]
  
  #check that makes sense
  stopifnot( sum(expressionset$visit_type == "Routine") == NCOL(expressionset))
  stopifnot(sum(duplicated(expressionset$patient_id)) == 0)
  return(expressionset)
}