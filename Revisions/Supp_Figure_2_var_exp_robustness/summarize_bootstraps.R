library(dplyr)
library(tidyr)
library(tibble)

IN_DIR <- "./snakemake_workflow/varpart_objects/"
files <- list.files(IN_DIR, recursive = TRUE) 

vp_objs <- lapply(files, function(f){
  readRDS(file.path(IN_DIR, f))
})

source("../../2021_07_08/scripts/util/Plotting/tbnk_featurename_replace.R")

datatypes <- unique(dirname(files))
names(datatypes) <- datatypes
datatype_vec <- dirname(files)
vp_objs <- lapply(datatypes, function(dtype){vp_objs[datatype_vec == dtype]})

vp_objs <- lapply(vp_objs, function(vp_list){
  lapply(vp_list, data.frame)
})

summ_list <- lapply(vp_objs, function(vp_list){
  vp_list %>% 
          lapply(rownames_to_column, "feature") %>%
          bind_rows(.id = "boot_number") %>%
          group_by(feature) %>%
          summarise(mean_vexp = mean(patient_id), var_vexp = var(patient_id),
          lower = quantile(patient_id, .025),
          upper = quantile(patient_id, .975)
          )
})

plot_list <- lapply(names(summ_list), function(nm){
  summ <- summ_list[[nm]]
  summ <- summ %>%
          mutate(feature = factor(feature, levels = feature[order(mean_vexp)]))

  if(nm == "array_mod"){
          levels(summ$feature) <- replace_mod_names_single_type(levels(summ$feature), 1, excelpath = "../../2021_07_08/Inputs/Reference/module_name_map.xlsx")
  }
  if(nm == "soma_mod"){
          levels(summ$feature) <- replace_mod_names_single_type(levels(summ$feature), 2, excelpath = "../../2021_07_08/Inputs/Reference/module_name_map.xlsx")
  }
  if(nm == "cbctbnk"){
          levels(summ$feature) <- replace_tbnk_names(levels(summ$feature), "../../2021_07_08/Inputs/Reference/tbnk_names.csv")
  }

  p <- ggplot(summ, aes(x = feature, y = mean_vexp)) +
          geom_linerange(aes(ymin = lower, ymax = upper)) + 
          geom_point(color = "red") +
          theme_classic() +
          ylab("Variance Explained by Subject")
          ggtitle(nm)

  if(nm %in% c("array", "soma")){
    p <- p + 
          theme_classic() +
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
            ylab("Variance Explained by Subject") +
          ggtitle(nm)
  }else{
    p <- p + 
          theme_classic() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
            ylab("Variance Explained by Subject") +
          ggtitle(nm)
  }
 
})

plot_list2 <- lapply(names(summ_list), function(nm){
  summ <- summ_list[[nm]]


  if(nm %in% c("array", "soma")){
    p <- ggplot(summ, aes(x = mean_vexp, y = var_vexp)) +
          geom_point() +
          theme_classic() +
          geom_smooth() +
          ggtitle(nm)
  }else{
    p <- ggplot(summ, aes(x = mean_vexp, y = var_vexp)) +
          geom_text(aes(label = feature)) +
          theme_classic() +
          geom_smooth() +
          ggtitle(nm)
  }
 
})

dir.create("plots")
pdf("plots/vexp_with_confintervals_all_datatypes.pdf")
for(p in plot_list){
  print(p)
}
for(p in plot_list2){
  print(p)
}
dev.off()


pdf("plots/vexp_with_confintervals_all_datatypes_smaller.pdf", height = 4, width = 4)
for(p in plot_list){
  print(p)
}
for(p in plot_list2){
  print(p)
}
dev.off()
