
files <- list.files(".", recursive =TRUE, pattern = ".R")

files <- files[grepl("Paper_Figures", files)]

libs_list <- list() 

for(i in seq_along(files)){
  f <- files[[i]]
  txt <- readLines(f)
  txt <- txt[grepl("library\\(", txt)]
  txt <- gsub("library\\(","", txt)
  txt <- gsub("\\)", "", txt)
  txt <- gsub(" ", "", txt)
  libs_list[[i]] <-txt
}

libs <- unique(unlist(libs_list))

