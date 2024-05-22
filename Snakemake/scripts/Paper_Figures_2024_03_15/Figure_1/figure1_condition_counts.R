library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(Biobase)

if(!exists("snakemake")){
  setwd("../../..")
  source("scripts/util/paper/parse_snakemake.R")
  parse_snakemake("figure_1c_condition_counts")
}

# Source utilities
source('scripts/util/paper/abbrev_cond.R')

# Set paths
META.DATA.IN.PATH = snakemake@input[["meta"]]#'Metadata/monogenic.de-identified.metadata.RData'
#META.DATA.IN.PATH = 'Metadata/monogenic.de-identified.metadata.RData'

COND.GROUP.IN.PATH <- snakemake@input[["cond_groups"]]#"Reference/condition_groups.csv"

FIGURE.1c.OUT.PATH = snakemake@output[[1]]#"Paper_1_Figures/Figure_1/1c_v2.pdf"
cond_group_dat <- read.csv(COND.GROUP.IN.PATH, stringsAsFactors = FALSE, check.names = FALSE)


# Figure 1c -- stacked barplots displaying counts for conditions, broken down by primary and repeated samples

## Load the monogenic metadata database
load(META.DATA.IN.PATH)

## Get the set of conditions represented in the Discovery and Validation sets (not QC sets)
#conditions = monogenic.all.assays[monogenic.all.assays$analysis_group %in% c('Discovery','Validation'), 'condition']
conditions = monogenic.all.assays[monogenic.all.assays$analysis_group %in% c('Discovery'), 'condition']

## Remove any unknown samples from these conditions
unknowns = c('Unknown', 'U_PI3K', 'U_Kastner', 'U_STAT1', 'U_CTLA4', 'U_Telomere')
conditions = setdiff(conditions, unknowns)

#train_test_count <- monogenic.all.assays %>%
#  filter(assay_type %in% c("Microarray", "Somalogic")) %>%
#  select(patient_id, visit_id, analysis_group, condition) %>%
#  filter(!grepl("CTRL", condition, ignore.case = T)) %>%
#  filter(!grepl("control", condition, ignore.case = T)) %>%
#  mutate(condition2 = condition) %>%
#  mutate(condition = replace(condition, !condition %in% conditions, "Other")) %>%
#  unique() %>% # Ensure there are no repeated visits
#  ##filter(analysis_group %in% 'Discovery') %>% # Filter to just training samples
#  select(patient_id, condition, analysis_group) %>%
#  unique() %>%
#  group_by(condition) %>%
#  summarise(ntrain = sum(analysis_group == "Discovery"), 
#  ntest= sum(analysis_group == "Validation")) %>% as.data.frame()
#
#sum(train_test_count$ntest)


#counts_other = monogenic.all.assays %>%
#  filter(assay_type %in% c("Microarray", "Somalogic")) %>%
#  select(patient_id, visit_id, analysis_group, condition) %>%
#  filter(!grepl("CTRL", condition, ignore.case = T)) %>%
#  filter(!grepl("control", condition, ignore.case = T)) %>%
#  mutate(condition2 = condition) %>%
#  mutate(condition = replace(condition, !condition %in% conditions, "Other")) %>%
#  unique() %>% # Ensure there are no repeated visits
#  ##filter(analysis_group %in% 'Discovery') %>% # Filter to just training samples
#  select(-analysis_group) %>%
#  group_by(condition2) %>%
#  summarise(primary.samples = length(unique(patient_id)), # Count the number of primary samples for each condition
#            other = sum(condition == "Other") != 0,
#            all.samples = length(unique(visit_id))) %>% # Count the number of total samples for each condition
#  mutate(repeat.samples = all.samples - primary.samples) %>% # Get the number of repeat samples for each condtion
#  ungroup() %>%
#  arrange(desc(all.samples)) %>% # Sort by the total number of samples
#  filter(other) %>%
#  select(-other)

## Create a matrix with counts of primary and repeat training set samples for each condition
counts_both_test_train = monogenic.all.assays %>%
  filter(assay_type %in% c("Microarray", "Somalogic")) %>%
  select(patient_id, visit_id, analysis_group, condition) %>%
  filter(!grepl("CTRL", condition, ignore.case = T)) %>%
  filter(!grepl("control", condition, ignore.case = T)) %>%
  #mutate(condition2 = condition) %>%
  #mutate(condition = replace(condition, !condition %in% conditions, "Other")) %>%
  filter(condition %in% c("Healthy", conditions)) %>%
  unique() %>% # Ensure there are no repeated visits
  ##filter(analysis_group %in% 'Discovery') %>% # Filter to just training samples
  group_by(condition, analysis_group) %>%
  summarise(primary.samples = length(unique(patient_id)), # Count the number of primary samples for each condition
            all.samples = length(unique(visit_id))) %>% # Count the number of total samples for each condition
  mutate(repeat.samples = all.samples - primary.samples) %>% # Get the number of repeat samples for each condtion
  ungroup() 

counts_train <- counts_both_test_train %>% 
        filter(analysis_group == "Discovery") %>%
        select(-c(analysis_group, all.samples))
counts_test <- counts_both_test_train %>% filter(analysis_group == "Validation") %>%
        select(condition, primary.samples) %>%
        rename(set.aside.samples =primary.samples)

counts <- left_join(counts_train, counts_test) %>%
        mutate(set.aside.samples = replace(set.aside.samples, is.na(set.aside.samples), 0))

counts <- counts %>%
  group_by(condition) %>%
  left_join(cond_group_dat) %>%
  arrange(desc(primary.samples + set.aside.samples + repeat.samples)) %>% # Sort by the total number of samples
  #select(-all.samples) %>%
  mutate(cond_abbrev = abbrev_cond(condition)) # Change the condition column to be a factor with levels in the same order

counts$cond_abbrev <- factor(counts$cond_abbrev, levels = rev(counts$cond_abbrev))


counts <- counts %>%
  gather(key = variable, value = value, -c(condition, cond_group, cond_abbrev)) %>%
  mutate(variable = as.character(variable)) %>%
  mutate(variable = replace(variable, variable == 'primary.samples', 'Primary Sample')) %>%
  mutate(variable = replace(variable, variable == 'repeat.samples', 'Additional Samples:\nCollected at different visits/timepoints')) %>%
  mutate(variable = replace(variable, variable == 'set.aside.samples', 'Set Aside Sample')) %>%
  mutate(variable = factor(variable, levels = rev(c('Primary Sample', "Set Aside Sample", 'Additional Samples:\nCollected at different visits/timepoints')))) %>%
  ungroup()


counts <- counts %>%
        mutate(cond_group = as.character(cond_group)) %>%
        mutate(cond_group = replace(cond_group, condition == "Healthy", "Healthy")) %>%
        mutate(cond_group = replace(cond_group, cond_group == "TERT.TERC", "Telo")) %>%
        mutate(cond_group = factor(cond_group, levels = c("Healthy", "AI", "Telo", "PID")))


## Create the barplots
#p = ggplot(counts, aes(x = cond_abbrev, y = value, fill = variable)) + 
p = ggplot(counts, aes(x = cond_abbrev, y = value, fill = cond_group, alpha = variable)) + 
        geom_bar(stat = 'identity') + 
  ylab('# Samples') + xlab('Condition') + labs(alpha = 'Sample Type', fill = "Group") + theme_bw() + 
  facet_grid(cond_group~1, scales = "free", space = "free") +
  theme(
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    strip.text = element_blank(),
    strip.background = element_blank(),
  ) + coord_flip() +
  scale_alpha_manual(values = c(.4, .8, 1))
   #scale_fill_manual(values = c('steelblue2','steelblue4'), guide = guide_legend(reverse = TRUE))

## Save the barplots
ggsave(FIGURE.1c.OUT.PATH, p, device = 'pdf', height = 7, width = 10)

