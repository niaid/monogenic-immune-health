library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character",
        bench_iteration = "numeric",
        scriptdir = "character",
        source = "function"
    )
)

#snakemake <- new("Snakemake")

parse_snakemake <- function(rule, snakefile_path = "Snakefile"){

  txt <- readLines(snakefile_path)
  
  i_start <- which(txt == paste0("rule ", rule, ":"))
  if(length(i_start)== 0){
    stop(paste0(rule, " is not a rule defined in the Snakefile"))
  }
  i_end <- length(txt)
  for(i in seq(i_start, length(txt))){
    if (txt[[i]] == ""){
       i_end = i -1
       break()
    }
    i <- i +1
  }
  state <- NULL
  
  params_vec <- vector("character")
  params_list <- list()
  for(i in seq(i_start +1, i_end)){
    txt_single <- txt[[i]]
    if(grepl("#", txt_single)){
      txt_single <- strsplit(txt_single, "#")[[1]][[1]]
    }
  
    txt_single <- gsub(" ", "", txt_single)
    txt_single <- gsub(",", "", txt_single)
    txt_single <- gsub("\\\"", "", txt_single)
    txt_single <- gsub("\\\'", "", txt_single)
    txt_single <- gsub("\t", "", txt_single)
  
    if(startsWith(txt_single, "#")){
      next()
    }
  
    if(endsWith(txt_single, ":")){
      state <- gsub(":", "", txt_single)
      params_list[[state]] <- vector("character")
      next()
    }
    if(grepl("=", txt_single)){
      splt <- strsplit(txt_single, split ="\\=")
      entry <- splt[[1]][[2]]
      names(entry) = splt[[1]][[1]]
    }else{
      entry = txt_single
    }
  
    params_list[[state]] <- c(params_list[[state]], entry)
  }

  snakemake_tmp <- new("Snakemake")
  for(p in names(params_list)){
    if(p == "script"){
      next()
    }
    slot(snakemake_tmp, p) <- as.list(params_list[[p]])
  }
  snakemake <<- snakemake_tmp
}

#to test function
#parse_snakemake("../Snakefile", "all")


#testing not as function -----
#path <- "../Snakefile"
#
#txt <- readLines(path)
#
#grep("rule", txt, value = T)
#
##rule <- "all"
##rule <- "run_gvi_permutations"
#rule <- "compile_supp_tables_excel"
#
#i_start <- which(txt == paste0("rule ", rule, ":"))
#i_end <- length(txt)
#for(i in seq(i_start, length(txt))){
#  if (txt[[i]] == ""){
#     i_end = i -1
#     break()
#  }
#  i <- i +1
#}
#state <- NULL
#
#params_vec <- vector("character")
#params_list <- list()
#for(i in seq(i_start +1, i_end)){
#  txt_single <- txt[[i]]
#  if(grepl("#", txt_single)){
#    txt_single <- strsplit(txt_single)[[1]][[1]]
#  }
#
#  txt_single <- gsub(" ", "", txt_single)
#  txt_single <- gsub(",", "", txt_single)
#  txt_single <- gsub("\\\"", "", txt_single)
#  txt_single <- gsub("\\\'", "", txt_single)
#  txt_single <- gsub("\t", "", txt_single)
#
#  if(startsWith(txt_single, "#")){
#    next()
#  }
#
#  if(endsWith(txt_single, ":")){
#    state <- gsub(":", "", txt_single)
#    params_list[[state]] <- vector("character")
#    next()
#  }
#  if(grepl("=", txt_single)){
#    splt <- strsplit(txt_single, split ="\\=")
#    entry <- splt[[1]][[2]]
#    names(entry) = splt[[1]][[1]]
#  }else{
#    entry = txt_single
#  }
#
#  params_list[[state]] <- c(params_list[[state]], entry)
#}
#
#snakemake <- new("Snakemake")
#for(p in names(params_list)){
#  if(p == "script"){
#    next()
#  }
#  slot(snakemake, p) <- as.list(params_list[[p]])
#}
