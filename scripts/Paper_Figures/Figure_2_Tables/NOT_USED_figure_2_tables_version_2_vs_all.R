# Load libraries
library(reshape2)
library(dplyr)
library(tibble)
library(Biobase)

# Set paths


setwd("/hpcdata/sg/sg_data/PROJECTS/Monogenic_Project/")
SOMALOGIC_MODULE_IN_PATH <- "Data/Somalogic/analysis_output/wgcna_results/modules.rds"
MICROARRAY_MODULE_IN_PATH <- "Data/Microarray/analysis_output/WGCNA/modules.rds"

SOMALOGIC_OUT_PATH <- "Paper_1_Figures/Figure_1_Tables/somalogic_module_members.csv"
MICROARRAY_OUT_PATH <- "Paper_1_Figures/Figure_1_Tables/microarray_module_members.csv"

soma_mods <- readRDS(SOMALOGIC_MODULE_IN_PATH)
array_mods <- readRDS(MICROARRAY_MODULE_IN_PATH)

source("scripts/util/Plotting/tbnk_featurename_replace.R")

soma_mods <- soma_mods %>% enframe(value = "module_color", name = "feature") %>%
        mutate(module_name = replace_mod_names_single_type(module_color, sheet = "PM")) %>%
        select(feature, module_name, module_color )

array_mods <- array_mods %>% enframe(value = "module_color", name = "feature") %>%
        mutate(module_name = replace_mod_names_single_type(module_color, sheet = "TM")) %>%
        select(feature, module_name, module_color)

## Limma DE statistics
RESULTS.IN.PATHS = list(
  `protein modules` ='Data/Somalogic/analysis_output/differential_expression/somalogic_modules_DE_results.rds',
  `protein features` = 'Data/Somalogic/analysis_output/differential_expression/somalogic_features_DE_results.rds',
  `gene modules` = 'Data/Microarray/analysis_output/differential_expression/microarray_modules_DE_results.rds',
  `gene features` = 'Data/Microarray/analysis_output/differential_expression/microarray_features_DE_results.rds',
  `cbcs and tbnks` = 'Data/TBNK/analysis_output/differential_expression/tbnks_features_DE_results.rds'
)

## Somalogic eset
PROTEIN.ESET.IN.PATH = 'Data/Somalogic/data_analysis_ready/analysis_ready_sample_level_training_somalogic.rds'

## Tables out paths
TABLE.OUT.PATHS = list(
  `protein modules` = 'Paper_1_Figures/Figure_2_Tables/protein_modules_DE_results_vs_all.txt',
  `protein features` = 'Paper_1_Figures/Figure_2_Tables/protein_features_DE_results_vs_all.txt',
  `gene modules` = 'Paper_1_Figures/Figure_2_Tables/gene_modules_DE_results_vs_all.txt',
  `gene features` = 'Paper_1_Figures/Figure_2_Tables/gene_features_DE_results_vs_all.txt',
  `cbcs and tbnks` = 'Paper_1_Figures/Figure_2_Tables/cbc_and_tbnks_DE_results_vs_all.txt'
)

# Load data
results = lapply(RESULTS.IN.PATHS, readRDS)
somalogic = readRDS(PROTEIN.ESET.IN.PATH)

# For each data type
results.long = lapply(names(results), function(data.type) {
  
  # Get the data for that data type
  result = results[[data.type]]
  
  # Decide whether the data type is a module or feature
  if(grepl('module', data.type)) {
    feature.header = 'module'
  } else {
    feature.header = 'feature'
  }
  
  # For each statistic
  stats = lapply(names(result$versus.all), function(stat) {
    
    # Convert statistics from wide to long
    df = melt(result$versus.all[[stat]])
    
    # Rename the columns
    colnames(df) = c(feature.header, 'condition', stat)
    
    # Add the data tpye
    df$data.type = data.type
    
    # Select the desired columns
    df = df %>% select(!!feature.header, condition, data.type, !!stat)
    return(df)
  })
  
  # Name the long form statistics after the data tpyes
  names(stats) = names(result$versus.all)
  
  # Combine the data frames to put each statistic into one data frame
  result.long = stats[[1]]
  for(stat in stats[2:length(stats)]) {
    result.long = result.long %>% right_join(stat, by = c(feature.header,'condition','data.type'))
  }
  
  return(result.long)
})

names(results.long) = names(results)

# Add the feature meta data from the proteins to the proteomic data frame as the feature names are not standardized like genes
protein.features = results.long$`protein features`
f = fData(somalogic)[protein.features$feature,]
protein.features = cbind(protein.features, f)
results.long$`protein features` = protein.features

results.long$`protein features` <- right_join(soma_mods, results.long$`protein features` )
results.long$`gene features` <- right_join(array_mods, results.long$`gene features` )

results.long$`protein modules` <- results.long$`protein modules` %>%
        rename(module_color = module) %>%
        mutate(module_name = replace_mod_names_single_type(as.character(module_color), sheet = "PM")) %>%
        select(module_name, module_color, everything())

results.long$`gene modules` <- results.long$`gene modules` %>%
        rename(module_color = module) %>%
        mutate(module_name = replace_mod_names_single_type(as.character(module_color), sheet = "TM")) %>%
        select(module_name, module_color, everything())

# Output tables to text files
for(data.type in names(results.long)) {
  result.long = results.long[[data.type]]
  file.path = TABLE.OUT.PATHS[[data.type]]
  write.table(result.long, file = file.path, sep = '\t', row.names = F, col.names = T)
}
