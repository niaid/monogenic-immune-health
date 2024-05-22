library(tidyverse)
library(reshape2)
library(cowplot)

if(!exists("snakemake")){
  setwd("../../..")
  source("scripts/util/paper/parse_snakemake.R")
  parse_snakemake(rule = "figure_1_tbnk_mod_vp")
}


MICROARRAY.MODULES.VP.IN.PATH = snakemake@input[["array_mod"]]#'Data/Microarray/analysis_output/variance_decomposition/microarray_modules_vp.RDS'
SOMALOGIC.MODULES.VP.IN.PATH = snakemake@input[["soma_mod"]]#'Data/Somalogic/analysis_output/variance_decomposition/somalogic_modules_vp.RDS'
TBNKS.VP.IN.PATH = snakemake@input[["tbnk"]]#'Data/TBNK/analysis_output/variance_decomposition/tbnk_features_vp.RDS'

FIG.OUT.PATH <- snakemake@output[[1]]#"Paper_1_Figures/Figure_1/module_tbnk_varpart.pdf"

#setwd("../../..")
#MICROARRAY.MODULES.VP.IN.PATH = 'Pipeline_out/Data/Microarray/analysis_output/variance_decomposition/microarray_modules_vp.RDS'
#SOMALOGIC.MODULES.VP.IN.PATH = 'Pipeline_out/Data/Somalogic/analysis_output/variance_decomposition/somalogic_modules_vp.RDS'
#TBNKS.VP.IN.PATH = 'Pipeline_out/Data/TBNK/analysis_output/variance_decomposition/tbnk_features_vp.RDS'
#
#FIG.OUT.PATH <- "Pipeline_out/Paper_1_Figures/Figure_1/module_tbnk_varpart.pdf"

source("scripts/util/Plotting/tbnk_featurename_replace.R")

# Figure 1f -- stacked bar plots of variance partitions for modules and tbnks

## Load the data
microarray.vp = readRDS(MICROARRAY.MODULES.VP.IN.PATH)
somalogic.vp = readRDS(SOMALOGIC.MODULES.VP.IN.PATH)
tbnks.vp = readRDS(TBNKS.VP.IN.PATH)

## Define the medications we wish to include in the summary variance for medication
medications = c('IgG.replacement', 'Anti.TNF', 'IFN.gamma', 'Immune.stimulators', 
                'Anti.IL1', 'Antifungal', 'Steroid', 'Anti.inflammatories',
                'Antibiotic', 'Immunosuppressant', 'Antibody')

## Instantiate a function to summarize each variance partition into a data frame
extract_results = function(results, medications) {
  
  ### Convert the variance parititon results to a data frame
  df = data.frame(results)
  
  ### Create the initial data frame, which displays the variation associated with patient and condition
  ### and which summarizes variance attributed to each medication type into a single score
  df = df %>% 
    tibble::rownames_to_column(var = 'module') %>%
    mutate(Medication = rowSums(as.matrix(df[, medications]))) %>%
    select(-!!medications) %>%
    select(-Residuals) %>%
    mutate(module = factor(module, levels = rev(sort(unique(module))))) %>%
    melt(id.vars = 'module', variable.name = 'Covariate') %>%
    mutate(Covariate = factor(Covariate, levels = rev(c('Patient', 'Condition', 'Medication'))))
  
  ### Get the amount of variation associated with patient for each feature
  df.patient = df %>% 
    filter(Covariate == 'Patient') %>%
    arrange(value)

  ### Arrange the modules to order by patient-associated variation
  df = df %>%
    mutate(module = module %>% as.character(module)) %>%
    mutate(module = factor(module, levels = unique(df.patient$module)))
  
  return(df)
}

## Instantiate a function to create the barplot from the extracted results
bar_plot = function(df, colors) {
  p = ggplot(df, aes(x = module, y = value, fill = Covariate)) +
    geom_bar(stat = 'identity', show.legend = TRUE) +
    theme_bw() + ylim(0,1) + coord_flip() + 
    scale_fill_manual(values = colors, guide = guide_legend(reverse = TRUE))
}

## Panel 1 -- TBNKs
### Manually rename the tbnk features names to make them clearer and more concside
df.tbnks = extract_results(tbnks.vp, medications)
#levels(df.tbnks$module) = levels(df.tbnks$module) %>% 
levels(df.tbnks$module) = levels(df.tbnks$module) %>% replace_tbnk_names()

df.tbnks <- df.tbnks %>%
        mutate(category = tbnk_groups(module, "new name"))


### Choose the tbnk plotting colors
#colors = c('seagreen1','seagreen3','seagreen4')
colors = c('aquamarine','seagreen3','seagreen4')

### Create the tbnk barplot
p.tbnks = bar_plot(df.tbnks, colors) + 
  #xlab('Major Peripheral \nImmune Parameters') +
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
df.somalogic = extract_results(somalogic.vp, medications)

### Choose somalogic plotting colors
colors = c('thistle1','plum1','violetred')

levels(df.somalogic$module) <- replace_mod_names_single_type(levels(df.somalogic$module), "PM")

### Create the somalogic barplot
p.somalogic = bar_plot(df.somalogic, colors) + 
  #xlab('Proteomic Modules') +
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
df.microarray = extract_results(microarray.vp, medications)
levels(df.microarray$module) <- replace_mod_names_single_type(levels(df.microarray$module), "TM")


### Choose the microarray colors
colors = c('steelblue2','royalblue','royalblue4')

### Create the microarray bar plot
p.microarray = bar_plot(df.microarray, colors) + 
  #xlab('Transcriptomic Modules') + 
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

ggsave(FIG.OUT.PATH, p, device = 'pdf', height = 12, width = 10)
