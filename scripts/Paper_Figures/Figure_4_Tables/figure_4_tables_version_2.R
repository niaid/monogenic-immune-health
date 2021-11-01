# Load library
library(dplyr)
library(Biobase)

# Set paths
if(!exists("snakemake")){
  setwd("../../..")
  source("scripts/util/paper/parse_snakemake.R")
  parse_snakemake(rule = "figure_4_tables")

}
GVIS.IN.PATHS = list(
  HEALTHY = snakemake@input[["gvis"]]#'Classification/results/healthy_rf_gvis_all.RDS',
)

PVALS.IN.PATHS = list(
  HEALTHY = snakemake@input[["pvals"]]#'Classification/results/healthy_rf_pvals_all.RDS',
)

SOMALOGIC.IN.PATH = snakemake@input[["soma_data"]]#'Data/Somalogic/data_analysis_ready/analysis_ready_sample_level_training_somalogic.rds'

TABLE.OUT.PATHS = list(
  HEALTHY = snakemake@output[["table"]]#'Paper_1_Figures/Figure_4_Tables/healthy_feature_gvi_table.txt',
)

source("scripts/util/Plotting/tbnk_featurename_replace.R")

# Create maps to rename the classifier groups and data types

classifier.groups = list(
  HEALTHY = 'Healthy versus all conditions'
)

classifier.names = list(
  'all.modules.plus.grey.with.tbnks' = 'CBCs, TBNKs, Gene Module Scores, Protein Module Scores, & Grey Module Proteins'
)

# Load the data
gviss = lapply(GVIS.IN.PATHS, readRDS)
pvalss = lapply(PVALS.IN.PATHS, readRDS)
eset = readRDS(SOMALOGIC.IN.PATH)

# Get the feature meta data from the eset and remove the Units column (it's somewhat confusing in the context of the other data types)
somamer.meta.data = fData(eset) %>% select(-Units)

# For each classification task (e.g. healthy versus all)
dfs = mapply(function(gvis, pvals, classifier.group) {
  # For each classifier (e.g. tbnks, protein modules, etc.)
  dfs = mapply(function(gvi, pval, classifier.name, classifier.number) {
    
    # Make sure the gvis and their corresponding pvalues are in the same order
    stopifnot(names(gvi) == names(pval))
    
    # Get the feature names associated with the gvis
    feature.names = names(gvi)
    
    # Initialize the data frame
    df = data.frame(`Feature Name` = feature.names, 
                    GVI = gvi, 
                    Pval = pval,
                    `Classifier Number` = classifier.number,
                    `Classifier Objective` = classifier.group,
                    `Data Types in Classifier` = classifier.name,
                    stringsAsFactors = F,
                    check.names = F)
    
    # Get the data type of the feature based on the header in the feature name
    df$`Feature Data Type` = df$`Feature Name` %>%
      replace(., grepl('^somalogic\\.grey\\.',.), 'Grey Module Proteins') %>%
      replace(., grepl('^tbnks\\.',.), 'CBCs/TBNKs') %>%
      replace(., grepl('^somalogic\\.modules\\.',.), 'Protein Module') %>%
      replace(., grepl('^microarray\\.modules\\.',.), 'Gene Module')
    
    # Manaully remove the feature name from the header
    df$`Feature Name` = df$`Feature Name` %>%
      gsub('^somalogic\\.grey\\.','',.) %>%
      gsub('^tbnks\\.','',.) 
    
    
    # Rearrang the data frame column orders
    df = df[, c('Classifier Number',
                'Classifier Objective',
                'Data Types in Classifier',
                'Feature Name',
                'Feature Data Type',
                'GVI',
                'Pval')]
    
    # Add the somamer metadata to help identify somamer
    df = cbind(df, somamer.meta.data[df$`Feature Name`,])
    
    return(df)
    
  }, gvis, pvals, classifier.names, 1:6, SIMPLIFY = F)
  df = Reduce(rbind, dfs)
}, gviss, pvalss, classifier.groups, SIMPLIFY = F)


dat <- dfs[[1]]
#Want to keep only classifier number 6 because that is the one that includes everything it seems
#table(dat$`Classifier Number`, dat$`Feature Data Type`)
#dup_feat <- dat[["Feature Name"]] %>% .[duplicated(.)] %>% unique()
#dat %>% filter(`Feature Name` %in% dup_feat)
#table(dat[["Data Types in Classifier"]])

dat <- dat %>% filter(`Classifier Number` == 6)

dat <- dat %>% 
        mutate(AdjP = p.adjust(Pval, method = "fdr"))
ix <- which(colnames(dat) == "Pval")
dat <- bind_cols(dat[, 1:ix], data.frame(AdjP = dat[["AdjP"]]), dat[, (ix + 1):(ncol(dat) -1)])

#lapply(dat, function(x){
#  if(length(unique(x)) > 5){
#          table(table(x))
#  }else{
#    table(x)
#  }
#})
#table(dat[["Feature Name"]])


dat <- dat %>% select(-1)

dat <- dat %>% mutate(`Feature Name` = replace_tbnk_names(`Feature Name`))

dat <- dat %>%
        mutate(#`Feature Name` = 
               `Feature Name`= 
               replace_mod_names_both(`Feature Name`, 
                                      proteome_prefix = "somalogic.modules.", 
                                      transcriptome_prefix = "microarray.modules.")) 

#dat %>% filter(grepl("odule", `Feature Data Type`, ignore.case = T), 
#  !grepl("grey", `Feature Data Type`, ignore.case = T)
#)

readr::write_csv(dat, TABLE.OUT.PATHS[[1]])
