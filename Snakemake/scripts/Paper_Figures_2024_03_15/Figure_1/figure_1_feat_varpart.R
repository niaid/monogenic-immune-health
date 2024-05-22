# Load libraries
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(reshape2)
library(Biobase)
library(cowplot)

# Source utilities
source('scripts/util/paper/abbrev_cond.R')

# Set paths
META.DATA.IN.PATH = snakemake@input[[1]]#'Metadata/monogenic.de-identified.metadata.RData'
MICROARRAY.FEATURES.VP.IN.PATH = snakemake@input[[2]]#'Data/Microarray/analysis_output/variance_decomposition/microarray_features_vp.RDS'
SOMALOGIC.FEATURES.VP.IN.PATH = snakemake@input[[3]]#'Data/Somalogic/analysis_output/variance_decomposition/somalogic_features_vp.RDS'

FIGURE.1e.OUT.PATH = snakemake@output[["e"]]#"Paper_1_Figures/Figure_1/1e.pdf"
FIGURE.1f.OUT.PATH = snakemake@output[["f"]]#"Paper_1_Figures/Figure_1/1f.pdf"

# Figure 1a -- cartoon

# Figure 1b -- cartoon, see box

# Figure 1d -- feature-feature correlation matrices organized by module for somalogic and microarray, with enrichment-based annotations
## See scripts/Paper_Figures/Figure_1/figure_1_addendum_version_(latest).R for the script used to make the matrix plots
## See Paper_1_Figures/Figure_1/1d_proteomic.png and Paper_1_Figures/Figure_1/1d_transcriptomic.png for the plots themselves.
## These plots were saved as png rather than pdf because they are very large and saving them in high resolution would be prohibitively
## slow and memory intensive.
## These tables are generated using scripts/Enrichments/analysis/write_enrichment_directories.R

# Figure 1e -- stable versus unstable parameter example
set.seed(80)

## Randomly generate stable trajectories for the 'feature' in two patients
x = .8 ## x is the inital point of the feature for the first patient
xs = x ## xs is the time course of the feature for the first patient
y = .2 ## y is the inital point of the feature for the second patient
ys = y ## ys is the time course of the feature for the second patient
for(i in 1:999) { ## To get from one time point to the next
  ## randomly perturb the feature, with an elastic force pulling it back to its original value
  xs[[i + 1]] = xs[[i]] + .01 * rnorm(1) + .05 * (x - xs[[i]]) 
  ys[[i + 1]] = ys[[i]] + .01 * rnorm(1) + .05 * (y - ys[[i]])
}

## Create a data frame combining the stable trajectories
df1 = data.frame(subject = factor(c(rep(1, 1000), rep(2, 1000))), time = c(1:1000, 1:1000), parameter = c(xs, ys),
                 group = 'stable parameter')

## Randomly generate unstable trajectories for the 'feature' in two patients
x = .52 ## x is the inital point of the feature for the first patient
xs = x ## xs is the time course of the feature for the first patient
y = .48 ## y is the inital point of the feature for the second patient
ys = y ## ys is the time course of the feature for the second patient
for(i in 1:999) { ## To get from one time point to the next
  ## randomly perturb the feature, with an elastic force pulling it back to its original value
  xs[[i + 1]] = xs[[i]] + .05 * rnorm(1) + .05 * (x - xs[[i]])
  ys[[i + 1]] = ys[[i]] + .05 * rnorm(1) + .05 * (y - ys[[i]])
}

## Create a data frame combining the stable trajectories
df2 = data.frame(subject = factor(c(rep(1, 1000), rep(2, 1000))), time = c(1:1000, 1:1000), parameter = c(xs, ys),
                 group = 'unstable parameter')

## Join the stable and unstable trajectory dataframes
df = rbind(df1, df2)

## Make a line plot for each trajectory, separating by stability
p = ggplot(df, aes(x = time, y = parameter, color = subject)) + 
  scale_color_manual(values = c('steelblue2','lightcoral')) + 
  geom_line(show.legend = FALSE) + theme_bw() + facet_wrap(~group, nrow = 2) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text = element_text(size = 20)) + ylab('parameter value')

## Save the plot
ggsave(FIGURE.1e.OUT.PATH, p, device = 'pdf', height = 6, width = 9)

# Figure 1e -- violin plots of feature level variance partitioning for microarray and somalogic data

## Load the variance partitions
microarray.vp = readRDS(MICROARRAY.FEATURES.VP.IN.PATH)
somalogic.vp = readRDS(SOMALOGIC.FEATURES.VP.IN.PATH)

## Define the medications we wish to include in the summary variance for medication
medications = c('IgG.replacement', 'Anti.TNF', 'IFN.gamma', 'Immune.stimulators', 
                'Anti.IL1', 'Antifungal', 'Steroid', 'Anti.inflammatories',
                'Antibiotic', 'Immunosuppressant', 'Antibody')

## Instantiate a function to summarize each variance partition into a data frame
summarize_vp = function(results, medications) {
  df = data.frame(results)
  df = df %>%
    tibble::rownames_to_column(var='module') %>%
    mutate(module = as.character(module)) %>%
    mutate(Medication = rowSums(as.matrix(df[, medications]))) %>%
    select(-!!medications) %>%
    melt(id.vars = 'module') %>%
    mutate(variable = factor(variable, levels = c('Patient','Condition','Medication','Residuals')))
  return(df)
}

## Instantiate a function to create violin plot
violin_plot = function(df, colors) {
  ggplot(df, aes(x = variable, y = value, fill = variable)) + theme_bw() + 
    geom_violin(scale = "width", position = position_dodge(.8), width = .7, show.legend = FALSE) +
    scale_fill_manual(values = colors) + ylab('Variance Explained') + xlab('Covariate')
}

## Get the summaried variance partitions
df.microarray = summarize_vp(microarray.vp, medications)
df.somalogic = summarize_vp(somalogic.vp, medications)

## Set colors for protein
colors = c('violetred','plum1','thistle1','grey')

## Make somalogic violin plot
p1 = violin_plot(df.somalogic, colors) + ggtitle('Proteomic Features') +
  theme(
    axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
    axis.title.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    title = element_text(size = 15))

## Set colors for microarray 
colors = c('royalblue4','royalblue','steelblue2','grey')

## Make microarray violin plot
p2 = violin_plot(df.microarray, colors) + ggtitle('Transcriptomic Features') +
  theme(
    axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
    axis.title.x = element_text(size = 15),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    title = element_text(size = 15))

## Combine the two plots together in a single grid
p = plot_grid(p1, p2, align = "h", ncol = 2, rel_widths = c(10,9))

## Save the results
ggsave(FIGURE.1f.OUT.PATH, p, device = 'pdf', height = 6, width = 9)

