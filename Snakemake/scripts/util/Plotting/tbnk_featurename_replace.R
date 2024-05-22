replace_tbnk_names <- function(vec, csvpath = "Inputs/Reference/tbnk_names.csv"){
  #tbnk_names_dat <- readr::read_csv("Inputs/Reference/tbnk_names.csv")
  tbnk_names_dat <- read.csv(csvpath, stringsAsFactors= FALSE,
                             check.names = FALSE)
  tbnk_names_dat <- tbnk_names_dat[tbnk_names_dat$old_name %in% vec, ]
  oldnames <- tbnk_names_dat[["old_name"]]
  newnames <- tbnk_names_dat[["new name"]]

  for(i in seq_along(oldnames)){
    oldname <- oldnames[[i]]
    newname <- newnames[[i]]
    vec[vec == oldname] <- newname
  }
  vec
}

tbnk_groups <- function(vec, which_col, csvpath = "Inputs/Reference/tbnk_names.csv"){
  #tbnk_names_dat <- readr::read_csv("Inputs/Reference/tbnk_names.csv")
  tbnk_names_dat <- read.csv(csvpath, stringsAsFactors= FALSE,
                             check.names = FALSE)
  tbnk_names_dat <- tbnk_names_dat[tbnk_names_dat[[which_col]] %in% vec, ]
  feat <- tbnk_names_dat[[which_col]]
  categ <- tbnk_names_dat[["category"]]
  #print(feat)
  #print(categ)

  categ[match(vec, feat)] 
}

replace_mod_names_single_type <- function(vec, sheet, excelpath = "Inputs/Reference/module_name_map.xlsx"){
  dat <- readxl::read_excel(excelpath, sheet = sheet)
  for(i in ncol(dat)){
    dat[[i]] <- as.character(dat[[i]])
  }
  dat <- dat[dat$old_name %in% vec, ]
  oldnames <- dat[["old_name"]]
  newnames <- dat[["new_name"]]

  for(i in seq_along(oldnames)){
    oldname <- oldnames[[i]]
    newname <- newnames[[i]]
    vec[vec == oldname] <- newname
  }
  vec
}

replace_mod_names_both <- function(vec, proteome_prefix = "protein.", transcriptome_prefix = "gene.", excelpath = "Inputs/Reference/module_name_map.xlsx"){
  library(dplyr)
  dat <- readxl::read_excel(excelpath, sheet = 3)
  for(i in ncol(dat)){
    dat[[i]] <- as.character(dat[[i]])
  }

  dat$old_name[dat$type == "TM"] <- paste0(transcriptome_prefix, dat$old_name[dat$type == "TM"])
  dat$old_name[dat$type == "PM"] <- paste0(proteome_prefix, dat$old_name[dat$type == "PM"])
  #dat <- dat[dat$old_name %in% vec, ]
  oldnames <- dat[["old_name"]]
  newnames <- dat[["new_name"]]

  for(i in seq_along(oldnames)){
    oldname <- oldnames[[i]]
    newname <- newnames[[i]]
    vec[vec == oldname] <- newname
  }
  vec
}


