abbrev_cond <- function(condition.vec){
  cond.abbrev <- as.character(condition.vec)
  cond.abbrev <- replace(cond.abbrev, cond.abbrev == "PAPA Syndrome", "PAPA")
  cond.abbrev <- replace(cond.abbrev, cond.abbrev == "CD14", "CARD14 DN")
  cond.abbrev <- replace(cond.abbrev, cond.abbrev == "STAT1 GOF", "STAT1 GOF")
  cond.abbrev <- replace(cond.abbrev, cond.abbrev == "Job", "STAT3 DN")
  cond.abbrev <- replace(cond.abbrev, cond.abbrev == "XCGD", "X-CGD")
  cond.abbrev <- replace(cond.abbrev, cond.abbrev == "47CGD", "p47-CGD")
  cond.abbrev <- replace(cond.abbrev, cond.abbrev == "PI3K", "APDS1")
  cond.abbrev <- replace(cond.abbrev, cond.abbrev == "Healthy", "Healthy")
  #cond.abbrev <- replace(cond.abbrev, cond.abbrev == "Muckle-Wells", "MW")
  
  cond.abbrev
}

abbrev_cond_2 <- function(condition.vec){
  cond.abbrev <- as.character(condition.vec)
  cond.abbrev <- replace(cond.abbrev, cond.abbrev == "PAPA Syndrome", "PAPA")
  cond.abbrev <- replace(cond.abbrev, cond.abbrev == "STAT1 GOF", "STAT1")
  cond.abbrev <- replace(cond.abbrev, cond.abbrev == "CD14", "CARD14 DN")
  cond.abbrev <- replace(cond.abbrev, cond.abbrev == "Job", "STAT3")
  cond.abbrev <- replace(cond.abbrev, cond.abbrev == "XCGD", "X-CGD")
  cond.abbrev <- replace(cond.abbrev, cond.abbrev == "47CGD", "p47-CGD")
  cond.abbrev <- replace(cond.abbrev, cond.abbrev == "PI3K", "APDS1")
  cond.abbrev <- replace(cond.abbrev, cond.abbrev == "Healthy", "Healthy")
  #cond.abbrev <- replace(cond.abbrev, cond.abbrev == "Muckle-Wells", "MW")
  
  cond.abbrev
}


group_cond <- function(vec, set_levels = c("Healthy", "AI", "Telo", "PID")){
  dat <- read.csv("Inputs/Reference/condition_groups.csv", stringsAsFactors= FALSE,
                             check.names = FALSE)
  dat <- dat[dat$condition %in% vec, ]
  oldnames <- dat[["condition"]]
  newnames <- dat[["cond_group"]]

  out <- newnames[match(vec, oldnames)]
  
  if(!is.null(set_levels)){
    out <- factor(out, levels = set_levels)
  }
}

#group_cond <- function(condition.vec){
#  condition.tab <- table(condition.vec)
#  large.condition <- names(condition.tab)[condition.tab > 10]
#  ai <- c('FMF','FCAS','Muckle-Wells','PAPA Syndrome','TRAPS','HIDS','DADA2')
#  telo <- c("TERT", "TERC")
#  cgd <- c("47CGD", "XCGD")
#  
#  condition.vec <- as.character(condition.vec)
#  condition2 <- as.character(condition.vec)
#  
#  condition2 <- replace(condition2, condition.vec %in% cgd, "CGD")
#  condition2 <- replace(condition2, condition.vec %in% ai, "Autoinflammatory")
#  condition2 <- replace(condition2, condition.vec %in% telo, "Telomere")
#  condition2 <- replace(condition2, !condition.vec %in% c(large.condition, ai, telo), "Other PID")
#  condition2 <- factor(condition2)
#  
#  condition2
#}
