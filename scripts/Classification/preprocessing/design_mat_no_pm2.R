
if(!exists("snakemake")){
  setwd("../../..")
  source("scripts/util/paper/parse_snakemake.R")
  parse_snakemake(rule = "design_mat_no_pm2")
  
}

IN.PATH <- snakemake@input[[1]]

OUT.PATH <- snakemake@output[[1]]

design_mat_list <- readRDS(IN.PATH)

#design_mat_list <- design_mat_list[c(4, 5, 6)]

#remove pm2 somalogic.modules.purple
for(nm in names(design_mat_list)){
  m <- design_mat_list[[nm]]
  m <- m[, setdiff(colnames(m), "somalogic.modules.purple")]
  design_mat_list[[nm]] <- m
}

saveRDS(design_mat_list, OUT.PATH)

