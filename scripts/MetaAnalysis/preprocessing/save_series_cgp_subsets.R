cgps <- readRDS("/Volumes/SG_DATA/PROJECTS/Monogenic_Project/Reference/cgps_clean.RDS")
series <- readRDS("/Volumes/SG_DATA/PROJECTS/Monogenic_Project/Reference/jamboree_series_matrix/series_matrix_list.rds")

out.dir <- "/Volumes/SG_DATA/PROJECTS/Monogenic_Project/Reference/series_subsets"
dir.create(out.dir)

write_single <- function(subjects, series.mat, path){
  x <- series.mat[series.mat$geo_accession %in% subjects, ]
  write.csv(x, file = path)
}

for(disease in names(cgps)){
  out.dir.sub <- file.path(out.dir, disease)
  dir.create(out.dir.sub)
  for(study.platform in names(cgps[[disease]])){
    out.dir.sub.sub <- file.path(out.dir.sub, study.platform)
    dir.create(out.dir.sub.sub)
    
    cases <- cgps[[disease]][[study.platform]]$cases
    controls <- cgps[[disease]][[study.platform]]$controls
    
    #print(file.path(out.dir.sub.sub, "cases.csv"))
    write_single(cases, series.mat = series[[study.platform]],
                 path = file.path(out.dir.sub.sub, "cases.csv"))

    write_single(controls, series.mat = series[[study.platform]],
                 path = file.path(out.dir.sub.sub, "controls.csv"))
  }
}
