# Load libraries
library(reshape2)
library(tidyverse)
library(Biobase)

# Set paths
if(exists("snakemake")){
  ENRICHMENTS.IN.PATHS = list(
    protein.gene.sets = snakemake@input[[1]],#'Data/Somalogic/analysis_output/enrichments/somalogic_module_gene_set_enrichments.RDS',
    protein.tissue.sets = snakemake@input[[2]],#'Data/Somalogic/analysis_output/enrichments/somalogic_module_tissue_set_enrichments.RDS',
    gene.gene.sets = snakemake@input[[3]]#'Data/Microarray/analysis_output/enrichments/microarray_module_gene_set_enrichments.RDS'
  )
  
  META.DATA.IN.PATH = snakemake@input[[4]]#'Metadata/monogenic.de-identified.metadata.RData'
  
  ENRICHMENT.TABLES.OUT.PATHS = list(
    protein.gene.sets = snakemake@output[[1]],#'Paper_1_Figures/Figure_1_Tables/somalogic_module_gene_set_enrichments_table.txt',
    protein.tissue.sets = snakemake@output[[2]],#'Paper_1_Figures/Figure_1_Tables/somalogic_module_tissues_set_enrichments_table.txt',
    gene.gene.sets = snakemake@output[[3]]#'Paper_1_Figures/Figure_1_Tables/microarray_module_gene_set_enrichments_table.txt'
  )
  
  DEMOGRAPHICS.TABLE.OUT.PATH = snakemake@output[[4]]#'Paper_1_Figures/Figure_1_Tables/demographics_table.txt'
}else{
  ENRICHMENTS.IN.PATHS = list(
    protein.gene.sets = 'Data/Somalogic/analysis_output/enrichments/somalogic_module_gene_set_enrichments.RDS',
    protein.tissue.sets = 'Data/Somalogic/analysis_output/enrichments/somalogic_module_tissue_set_enrichments.RDS',
    gene.gene.sets = 'Data/Microarray/analysis_output/enrichments/microarray_module_gene_set_enrichments.RDS'
  )
  
  META.DATA.IN.PATH = 'Metadata/monogenic.de-identified.metadata.RData'
  
  ENRICHMENT.TABLES.OUT.PATHS = list(
    protein.gene.sets = 'Paper_1_Figures/Figure_1_Tables/somalogic_module_gene_set_enrichments_table.txt',
    protein.tissue.sets = 'Paper_1_Figures/Figure_1_Tables/somalogic_module_tissues_set_enrichments_table.txt',
    gene.gene.sets = 'Paper_1_Figures/Figure_1_Tables/microarray_module_gene_set_enrichments_table.txt'
  )
  
  DEMOGRAPHICS.TABLE.OUT.PATH = 'Paper_1_Figures/Figure_1_Tables/demographics_table.txt'

  setwd("/hpcdata/sg/sg_data/PROJECTS/Monogenic_Project")
}

source("scripts/util/Plotting/tbnk_featurename_replace.R")

# Create tables for the enrichment modules
## Load data
enrichments.list = lapply(ENRICHMENTS.IN.PATHS, readRDS)

## For each type of enrichment
for(enrichment.name in names(enrichments.list)) {
  
  ## Get the enrichments corresponding to that enrichment type
  enrichments = enrichments.list[[enrichment.name]]
  
  ## Get the table output path for that enrichment type
  enrichment.out.path = ENRICHMENT.TABLES.OUT.PATHS[[enrichment.name]]
  
  ## For each module in those enrichments
  enrichments = lapply(names(enrichments), function(module) {
    
    ## Get the enrichments corresponding to that module
    enrichment = enrichments[[module]]
    
    ## Set the name of the enrichments to the gene set names
    enrichment$Set.Name = rownames(enrichment)
    
    ## Set the module to the module color
    enrichment$module = module
    
    ## Rearrang the appended columns to be the first two columns
    n = ncol(enrichment)
    enrichment = enrichment[,c(n, n-1, 1:(n-2))]
  })
  
  ## Combine all the modules' enrichment data frames
  enrichments = Reduce(rbind, enrichments)

  enrichments <- enrichments %>% rename(module_color = module)

  if(startsWith(enrichment.name, "protein")){
    enrichments <- enrichments %>% 
            mutate(module_name = replace_mod_names_single_type(module_color, sheet = "PM"))
  }else if(startsWith(enrichment.name, "gene")){
    enrichments <- enrichments %>% 
            mutate(module_name = replace_mod_names_single_type(module_color, sheet = "TM"))
  }

  enrichments <- enrichments[, c("module_name", setdiff(colnames(enrichments), "module_name"))]

  if(enrichment.name != "protein.tissue.sets"){
    enrichments <- enrichments %>%
            group_by(module_name) %>%
            mutate(rank = rank(across.Adjusted.Pvalue)) %>%
            filter(rank < 101) %>%
            arrange(across.Adjusted.Pvalue) %>%
            select(-rank)
  }

  enrichments <- enrichments %>%
        select(-Adjusted.Pvalue) %>%
        rename(Adjusted.Pvalue = across.Adjusted.Pvalue) %>%
        mutate(tmp = gsub("grey", "PM9999", module_name)) %>%
        arrange(tmp) %>%
        select(-tmp)

  ## Output the data frames as a table 
  write.table(enrichments, file = enrichment.out.path, row.names = FALSE, sep = '\t', quote = FALSE)
}

# Create tables for demographics
## Load Data
load(META.DATA.IN.PATH)

## ## Get metadata for patients in training set
## Get ages for patients in training set
df = monogenic.all.assays %>%
  select(visit_id, patient_id, patient_age_at_time_of_blood_draw, condition, analysis_group, gender, race, ethnicity) %>%
  mutate(patient_id = paste0('P', patient_id)) %>%
  mutate(age = patient_age_at_time_of_blood_draw) %>%
  mutate(group = ifelse(condition == "Healthy", "Control", "Case")) %>%
  filter(analysis_group == 'Discovery') %>%
  select(-patient_age_at_time_of_blood_draw, -analysis_group)

## Average ages over various visits, first by averaging ages from samples within a visit (should all be almost exactly the same),
## and then averaging ages from samples across visits
df = df %>%
  group_by(visit_id) %>%
  summarise(age = mean(age),
            patient_id = unique(patient_id),
            condition = unique(condition), 
            group = unique(group), 
            gender = unique(gender),
            ethnicity = unique(ethnicity),
            race = unique(race)) %>%
  ungroup() %>%
  group_by(patient_id) %>%
  summarise(age = mean(age),
            condition = unique(condition), 
            group = unique(group), 
            gender = unique(gender),
            ethnicity = unique(ethnicity),
            race = unique(race),
            visits = length(visit_id)) %>%
  ungroup()

## Get the statistics associated with each condition 
## (NOTE: I did not include the IQR/variance for age here so that patient-specific ages
## could not be estimated for conditions with 2 samples)
df = df %>%
  group_by(condition) %>%
  summarise(`number of patients` = length(patient_id),
            `number of samples` = sum(visits),
            `age median` = median(age),
            `percent Female (Sex)` = mean(gender == 'F'),
            `percent Male (Sex)` = mean(gender == 'M'),
            `percent Asian (Race)` = mean(race == 'Asian'),
            `percent Black / African American (Race)` = mean(race == "Black/African Amer"),
            `percent Hawaiian/Pacific Islander (Race)` = mean(race == "Hawaiian/Pac. Island"),
            `percent Unknown (Race)` = mean(race == "Unknown"),
            `percent Multi-Racial (Race)` = mean(race == "Multiple Race"),
            `percent White / Caucasian (Race)` = mean(race == "White"),
            `percent Hispanic or Latino (Ethnicity)` = mean(ethnicity == "Hispanic or Latino"),
            `percent Not Hispanic or Latino (Ethnicity)` = mean(ethnicity == "Not Hispanic or Latino"),
            `percent Unknown (Ethnicity)` = mean(ethnicity == "Unknown")) %>%
  mutate(condition = condition %>% factor(., sort(unique(.))) %>% relevel(., 'Healthy')) %>%
  arrange(condition) %>%
  mutate(`age median` = round(`age median`)) %>%
  as.data.frame()

## Print to table
write.table(df, DEMOGRAPHICS.TABLE.OUT.PATH, sep = '\t', row.names = F, col.names = T)
