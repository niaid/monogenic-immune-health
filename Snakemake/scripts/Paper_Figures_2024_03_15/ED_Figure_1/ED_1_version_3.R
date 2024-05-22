# Load libraries
library(Biobase)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(cowplot)

# Set paths
## Monogenic metadata
#if(exists("snakemake")){
METADATA.IN.PATH = snakemake@input[[1]]#'Metadata/monogenic.de-identified.metadata.RData'
## Microarray variance partition results with condition and medication covariates
MICROARRAY.FEATURE.VP.IN.PATH = snakemake@input[[2]]#'Data/Microarray/analysis_output/variance_decomposition/microarray_features_vp.RDS'
## Somalogic variance partition results with condition and medication covariates
SOMALOGIC.FEATURE.VP.IN.PATH = snakemake@input[[3]]#'Data/Somalogic/analysis_output/variance_decomposition/somalogic_features_vp.RDS'

## Microarray modules variance partition results with just patient as a covariate
MICROARRAY.MODULES.SIMPLE.VP.IN.PATH = snakemake@input[[4]]#'Data/Microarray/analysis_output/stability/microarray_modules_standard_vp.rds'
## Somalogic modules variance partition results with just patient as a covariate
SOMALOGIC.MODULES.SIMPLE.VP.IN.PATH = snakemake@input[[5]]#'Data/Somalogic/analysis_output/stability/somalogic_modules_standard_vp.rds'
## TBNKs modules variance partition results with just patient as a covariate
TBNKS.SIMPLE.VP.IN.PATH = snakemake@input[[6]]#'Data/TBNK/analysis_output/stability/tbnks_features_standard_vp.rds'
## Microarray features variance partition results with just patient as a covariate
MICROARRAY.FEATURES.SIMPLE.VP.IN.PATH = snakemake@input[[7]]#'Data/Microarray/analysis_output/stability/microarray_features_standard_vp.rds'
## Somalogic features variance partition results with just patient as a covariate
SOMALOGIC.FEATURES.SIMPLE.VP.IN.PATH = snakemake@input[[8]]#'Data/Somalogic/analysis_output/stability/somalogic_features_standard_vp.rds'

SUPPLEMENTAL.FIGURE.1a.OUT.PATH = snakemake@output[[1]]#'Paper_1_Figures/Supplemental_Figure_1/S1a.pdf'
SUPPLEMENTAL.FIGURE.1b.OUT.PATH = snakemake@output[[2]]#'Paper_1_Figures/Supplemental_Figure_1/S1b.pdf'
SUPPLEMENTAL.FIGURE.1c.OUT.PATH = snakemake@output[[3]]#'Paper_1_Figures/Supplemental_Figure_1/S1c.pdf'
SUPPLEMENTAL.FIGURE.1d.OUT.PATH = snakemake@output[[4]]#'Paper_1_Figures/Supplemental_Figure_1/S1d.pdf'
SUPPLEMENTAL.FIGURE.1e.OUT.PATH = snakemake@output[[5]]#'Paper_1_Figures/Supplemental_Figure_1/S1e.pdf'
SUPPLEMENTAL.FIGURE.1f.OUT.PATH = snakemake@output[[6]]#'Paper_1_Figures/Supplemental_Figure_1/S1f.pdf'
#}else{
#METADATA.IN.PATH = 'Metadata/monogenic.de-identified.metadata.RData'
### Microarray variance partition results with condition and medication covariates
#MICROARRAY.FEATURE.VP.IN.PATH = 'Data/Microarray/analysis_output/variance_decomposition/microarray_features_vp.RDS'
### Somalogic variance partition results with condition and medication covariates
#SOMALOGIC.FEATURE.VP.IN.PATH = 'Data/Somalogic/analysis_output/variance_decomposition/somalogic_features_vp.RDS'
#
### Microarray modules variance partition results with just patient as a covariate
#MICROARRAY.MODULES.SIMPLE.VP.IN.PATH = 'Data/Microarray/analysis_output/stability/microarray_modules_standard_vp.rds'
### Somalogic modules variance partition results with just patient as a covariate
#SOMALOGIC.MODULES.SIMPLE.VP.IN.PATH = 'Data/Somalogic/analysis_output/stability/somalogic_modules_standard_vp.rds'
### TBNKs modules variance partition results with just patient as a covariate
#TBNKS.SIMPLE.VP.IN.PATH = 'Data/TBNK/analysis_output/stability/tbnks_features_standard_vp.rds'
### Microarray features variance partition results with just patient as a covariate
#MICROARRAY.FEATURES.SIMPLE.VP.IN.PATH = 'Data/Microarray/analysis_output/stability/microarray_features_standard_vp.rds'
### Somalogic features variance partition results with just patient as a covariate
#SOMALOGIC.FEATURES.SIMPLE.VP.IN.PATH = 'Data/Somalogic/analysis_output/stability/somalogic_features_standard_vp.rds'
#
#SUPPLEMENTAL.FIGURE.1a.OUT.PATH = 'Paper_1_Figures/Supplemental_Figure_1/S1a.pdf'
#SUPPLEMENTAL.FIGURE.1b.OUT.PATH = 'Paper_1_Figures/Supplemental_Figure_1/S1b.pdf'
#SUPPLEMENTAL.FIGURE.1c.OUT.PATH = 'Paper_1_Figures/Supplemental_Figure_1/S1c.pdf'
#SUPPLEMENTAL.FIGURE.1d.OUT.PATH = 'Paper_1_Figures/Supplemental_Figure_1/S1d.pdf'
#SUPPLEMENTAL.FIGURE.1e.OUT.PATH = 'Paper_1_Figures/Supplemental_Figure_1/S1e.pdf'
#SUPPLEMENTAL.FIGURE.1f.OUT.PATH = 'Paper_1_Figures/Supplemental_Figure_1/S1f.pdf'
#}
# Source utilities
source('scripts/util/Plotting/colors.R')
source('scripts/util/Processing/averageRepeatSamples.R')
source('scripts/util/paper/abbrev_cond.R')
source("scripts/util/Plotting/tbnk_featurename_replace.R")


# Supplemental Figure 1a -- density plot of ages in case and control with lines for medians

## Load metadata
load(METADATA.IN.PATH)

## Get ages for patients in training set
df = monogenic.all.assays %>%
  select(visit_id, patient_id, patient_age_at_time_of_blood_draw, condition, analysis_group) %>%
  mutate(patient_id = paste0('P', patient_id)) %>%
  mutate(age = patient_age_at_time_of_blood_draw) %>%
  mutate(group = ifelse(condition == "Healthy", "Healthy", "Disease")) %>%
  mutate(group = factor(group, levels = c("Healthy", "Disease"))) %>%
  filter(analysis_group == 'Discovery') %>%
  select(-patient_age_at_time_of_blood_draw, -analysis_group)

## Average ages over various visits, first by averaging ages from samples within a visit (should all be almost exactly the same),
## and then averaging ages from samples across visits
df = df %>%
  group_by(visit_id) %>%
  summarise(age = mean(age), patient_id = unique(patient_id), condition = unique(condition), group = unique(group)) %>%
  ungroup() %>%
  group_by(patient_id) %>%
  summarise(age = mean(age), condition = unique(condition), group = unique(group)) %>%
  ungroup()

## Find median age for case and control
df.medians = df %>%
  group_by(group) %>%
  summarise(`Group Median Age` = median(age))

## Make density plot
p = ggplot(df, aes(x = age, fill = group)) + geom_density(alpha = .5) +
  geom_vline(aes(xintercept = `Group Median Age`, color = group), data = df.medians) +
  xlab('Age') + ylab('Density') + ggtitle('Age Distributions') +
  theme_bw() +
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 15))
ggsave(SUPPLEMENTAL.FIGURE.1a.OUT.PATH, p, device = 'pdf', height = 4, width = 6)

# Supplemental Figure 1b -- barplots of ages by condition

## Load metadata
load(METADATA.IN.PATH)

## Get metadata for patients in training set
## Get ages for patients in training set
df = monogenic.all.assays %>%
  select(visit_id, patient_id, patient_age_at_time_of_blood_draw, condition, analysis_group) %>%
  mutate(patient_id = paste0('P', patient_id)) %>%
  mutate(age = patient_age_at_time_of_blood_draw) %>%
  mutate(group = ifelse(condition == "Healthy", "Healthy", "Disease")) %>%
  mutate(group = factor(group, levels = c("Healthy", "Disease"))) %>%
  filter(analysis_group == 'Discovery') %>%
  select(-patient_age_at_time_of_blood_draw, -analysis_group)

## Average ages over various visits, first by averaging ages from samples within a visit (should all be almost exactly the same),
## and then averaging ages from samples across visits
df = df %>%
  group_by(visit_id) %>%
  summarise(age = mean(age), patient_id = unique(patient_id), condition = unique(condition), group = unique(group)) %>%
  ungroup() %>%
  group_by(patient_id) %>%
  summarise(age = mean(age), condition = unique(condition), group = unique(group)) %>%
  ungroup()

## Get the median age associated with each condition
median.ages = df %>%
  group_by(condition) %>%
  summarise(median.age = median(age)) %>%
  ungroup()

## Append the median age of the condition to the data frame, and order by median age
df = df %>%
  right_join(median.ages, by = 'condition') %>%
  arrange(median.age) %>%
  mutate(condition = condition %>% as.character %>% abbrev_cond) %>%
  mutate(condition = factor(condition, levels = unique(condition)))

## Create boxplots
p = ggplot(df, aes(x = condition, y = age, fill = group)) + 
  geom_boxplot(outlier.shape = NA) + 
  xlab('Condition') + ylab('Age') + ggtitle('Condition-Specific Age Distributions') +
  theme_bw() +
  theme(axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 15, vjust = .4),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 15))
ggsave(SUPPLEMENTAL.FIGURE.1b.OUT.PATH, p, device = 'pdf', height = 4, width = 6)

# Supplemental Figure 1c -- Gender split stacked barplot

## Load metadata
load(METADATA.IN.PATH)

## Get metadata for patients in training set
df = monogenic.all.assays %>%
  select(patient_id, gender, condition, analysis_group) %>%
  unique() %>%
  mutate(patient_id = paste0('P', patient_id)) %>%
  filter(analysis_group == 'Discovery')

## Get the total number of subjects of each gender and each condition
df = df %>%
  group_by(condition, gender) %>%
  summarise(gender.total = length(patient_id)) %>%
  ungroup()

## Get the total number of subjects of each condition
df.total = df %>%
  group_by(condition) %>%
  summarise(total = sum(gender.total)) %>%
  ungroup()

## Get the percent of subjects from each gender within a condition and sort by that fraction
df = df %>%
  right_join(df.total, by = 'condition') %>%
  mutate(percent = gender.total / total) %>%
  select(-gender.total, -total) %>%
  mutate(group = ifelse(condition == 'Healthy', 'Control', 'Case')) %>%
  mutate(percent.female = ifelse(gender == 'F', percent, 1 - percent)) %>%
  arrange(desc(group), desc(percent.female)) %>%
  mutate(condition = condition %>% as.character %>% abbrev_cond) %>%
  mutate(condition = factor(condition, levels = unique(condition)))

## Create the barplots
p = ggplot(df, aes(x = condition, y = percent, fill = gender)) + geom_bar(stat = 'identity') + 
  theme_bw() +
  xlab('Condition') + ylab('Percent') + ggtitle('Gender Split by Condition') +
  scale_fill_viridis_d() +
  theme(axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 15, vjust = .4),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 15))
ggsave(SUPPLEMENTAL.FIGURE.1c.OUT.PATH, p, device = 'pdf', height = 4, width = 6)

# Supplemental Figure 1d -- violin plots of medication-specific effects

## Extract protein and gene variance parititon results (with condition and medication covariates included)
microarray.vp = readRDS(MICROARRAY.FEATURE.VP.IN.PATH)
somalogic.vp = readRDS(SOMALOGIC.FEATURE.VP.IN.PATH)

## Extract the medication names
medications = setdiff(colnames(microarray.vp), c('Patient','Condition', 'Residuals'))

## Insantiate a function to summarize the variance partition into a dataframe
summarize_vp = function(results) {
  df = data.frame(results)
  df = melt(df)
  df$variable = factor(df$variable, levels = c('Patient','Condition', medications, 'Residuals'))
  return(df)
}

## Insantiate a function to make the violin plot from the extracted data frame
violin_plot = function(df, colors) {
  ggplot(df, aes(x = variable, y = value, fill = variable)) + theme_bw() + 
    geom_violin(scale = "width", position = position_dodge(.8), width = .7, show.legend = FALSE) +
    scale_fill_manual(values = colors) + ylab('Variance Explained') + xlab('Covariate')
}

## Extract the variance parititon results for the genes and proteins into separate data frames
df.microarray = summarize_vp(microarray.vp)
df.somalogic = summarize_vp(somalogic.vp)

## Make the plot for the proteins
colors = c('violetred','plum1', rep('thistle1', length(medications)),'grey')
p1 = violin_plot(df.somalogic, colors) + ggtitle('Proteomic Features') +
  theme(
    axis.text.x = element_text(size = 10, angle = 30, hjust = 1),
    axis.title.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10))

## Make the plot for the genes
colors = c('royalblue4','royalblue', rep('lightblue', length(medications)),'grey')
p2 = violin_plot(df.microarray, colors) + ggtitle('Transcriptomic Features') +
  theme(
    axis.text.x = element_text(size = 10, angle = 30, hjust = 1),
    axis.title.x = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank())

p = plot_grid(p1, p2, align = "h", ncol = 2, rel_widths = c(10,9))
ggsave(SUPPLEMENTAL.FIGURE.1d.OUT.PATH, p, device = 'pdf', height = 3, width = 7)

# Supplemental Figure 1e - Simple variance partition module and tbnk effects

## Load the data
microarray.vp = readRDS(MICROARRAY.MODULES.SIMPLE.VP.IN.PATH)
somalogic.vp = readRDS(SOMALOGIC.MODULES.SIMPLE.VP.IN.PATH)
tbnks.vp = readRDS(TBNKS.SIMPLE.VP.IN.PATH)

## Instantiate a function to summarize each variance partition into a data frame
extract_results = function(results) {
  df = data.frame(results)
  df = df %>% 
    select(-Residuals) %>%
    tibble::rownames_to_column(var = 'module') %>%
    arrange(Patient) %>%
    mutate(module = factor(module, levels = unique(module)))
  return(df)
}

## Instantiate a function to create the barplot from the extracted results
bar_plot = function(df, color) {
  p = ggplot(df, aes(x = module, y = Patient)) +
    geom_bar(stat = 'identity', fill = color, show.legend = TRUE) +
    theme_bw() + ylim(0,1) + coord_flip() + 
    geom_hline(yintercept = .5, linetype = 'dashed', color = 'black') +
    ylab('Percent variation explained by Patient')
}

## Panel 1 -- TBNKs
### Rename the tbnk features names to make them clearer and more concside
df.tbnks = extract_results(tbnks.vp)
levels(df.tbnks$module) = levels(df.tbnks$module) %>% replace_tbnk_names()

df.tbnks <- df.tbnks %>%
        mutate(category = tbnk_groups(module, "new name"))

### Create the tbnk barplot
p.tbnks = bar_plot(df.tbnks, color = 'seagreen3') + 
  xlab('') +
  facet_grid(category~1, space = "free", scales = "free_y") +
  ylab('Variance Explained') + 
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    strip.text.x = element_blank(),
    strip.background.x = element_blank()
  )

## Panel 2
df.somalogic = extract_results(somalogic.vp)

levels(df.somalogic$module) <- replace_mod_names_single_type(levels(df.somalogic$module), sheet = "PM")

### Create the somalogic barplot
p.somalogic = bar_plot(df.somalogic, color = 'violetred') + 
  xlab('') + 
  facet_grid("PM" ~ 1) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    strip.text.x = element_blank(),
    strip.background.x = element_blank()
  )

## Panel 3
df.microarray = extract_results(microarray.vp)
levels(df.microarray$module) <- replace_mod_names_single_type(levels(df.microarray$module), sheet = "TM")

### Create the microarray bar plot
p.microarray = bar_plot(df.microarray, color = 'royalblue4') + 
  xlab('') + 
  ylab('') + 
  facet_grid("TM" ~ 1) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    strip.text.x = element_blank(),
    strip.background.x = element_blank()
  )

## Put the panels together
p = plot_grid(p.microarray, p.somalogic, p.tbnks, 
              align = "hv", 
              axis = "tblr",
              nrow = 3, 
              rel_heights = c(nrow(df.microarray), nrow(df.somalogic) + 2, nrow(df.tbnks) + 7))
ggsave(SUPPLEMENTAL.FIGURE.1e.OUT.PATH, p, device = 'pdf', height = 12, width = 7)


# Supplemental Figure 1f - Feature percentiles for protein and gene
microarray.vp = readRDS(MICROARRAY.FEATURES.SIMPLE.VP.IN.PATH)
somalogic.vp = readRDS(SOMALOGIC.FEATURES.SIMPLE.VP.IN.PATH)

## Create the dataframes listing what percentile each feature is for variance explained by patient
## and what the corresponding variance explained is
df.microarray = microarray.vp %>%
  data.frame() %>%
  tibble::rownames_to_column(var = 'feature') %>%
  select(-Residuals) %>%
  mutate(data.type = 'WB Transcriptome') %>%
  arrange(Patient) %>%
  tibble::rowid_to_column(var = 'percentile') %>%
  mutate(percentile = percentile / nrow(microarray.vp))

df.somalogic = somalogic.vp %>%
  data.frame() %>%
  tibble::rownames_to_column(var = 'feature') %>%
  select(-Residuals) %>%
  mutate(data.type = 'Serum Proteins') %>%
  arrange(Patient) %>%
  tibble::rowid_to_column(var = 'percentile') %>%
  mutate(percentile = percentile / nrow(somalogic.vp))

## Put together the gene and protein data frames
df = rbind(df.microarray, df.somalogic)

## Create the percentile plots
p = ggplot(df, aes(x = percentile, y = Patient)) + 
  geom_bar(stat = 'identity', color = 'cyan', fill = 'cyan') + 
  facet_wrap(~ data.type, scales = 'free_x') +
  theme_bw() +
  labs(x = 'Percentile', y = 'Variance Explained')

ggsave(SUPPLEMENTAL.FIGURE.1f.OUT.PATH, p, device = 'pdf', height = 2.5, width = 5)
