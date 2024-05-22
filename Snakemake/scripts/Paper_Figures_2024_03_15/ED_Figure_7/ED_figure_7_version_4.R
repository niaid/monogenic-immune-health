# Load libraries
library(ggplot2)
library(Biobase)
library(ggpubr)
library(MetaIntegrator)
library(limma)
library(dplyr)


# Set paths
if(exists("snakemake")){
  ## The effect sizes associated with each feature in the meta analysis
  META.ANALYSIS.Z.SCORE.IN.PATH = snakemake@input[[1]]#'Reference/jamboree/analysis_output/results/jamboree_z_score_results.RDS'
  ## The comparison group pairs from the jamboree
  CGPS.IN.PATH = snakemake@input[[2]]#'Reference/jamboree/data_analysis_ready/cgps_clean.RDS'
  ## The transcriptional surrogate signature gene set enrichment results associated with each feature in the meta-analysis
  META.ANALYSIS.ENRICHMENTS.IN.PATH = snakemake@input[[3]]#'Reference/jamboree/analysis_output/results/jamboree_enrichment_results.RDS'
  ## The baltimore aging cohort eset
  ESET.IN.PATH = snakemake@input[[4]]#'Reference/ferrucci/processed/aging_eset.RDS'
  ## The plasma somalogic surrogate signature for the healthy index
  SIGNATURE.IN.PATH = snakemake@input[[5]]#'Classification/proteomic_surrogates/healthy.index.plasma.surrogate.somaId.RDS'
  ## The Ferrucci significance table
  TABLE.IN.PATH = snakemake@input[[6]]#'Reference/ferrucci/raw/acel12799-sup-0004-TableS3.txt'
  
  SUPPLEMENTAL.FIGURE.5a.OUT.PATH = snakemake@output[[1]]#'Paper_1_Figures/Supplemental_Figure_5/S5a.pdf'
  SUPPLEMENTAL.FIGURE.5b.OUT.PATH = snakemake@output[[2]]#'Paper_1_Figures/Supplemental_Figure_5/S5b.pdf'
  SUPPLEMENTAL.FIGURE.5c.OUT.PATH = snakemake@output[[3]]#'Paper_1_Figures/Supplemental_Figure_5/S5c.pdf'
}else{
  ## The effect sizes associated with each feature in the meta analysis
  META.ANALYSIS.Z.SCORE.IN.PATH = 'Reference/jamboree/analysis_output/results/jamboree_z_score_results.RDS'
  ## The comparison group pairs from the jamboree
  CGPS.IN.PATH = 'Reference/jamboree/data_analysis_ready/cgps_clean.RDS'
  ## The transcriptional surrogate signature gene set enrichment results associated with each feature in the meta-analysis
  META.ANALYSIS.ENRICHMENTS.IN.PATH = 'Reference/jamboree/analysis_output/results/jamboree_enrichment_results.RDS'
  ## The baltimore aging cohort eset
  ESET.IN.PATH = 'Reference/ferrucci/processed/aging_eset.RDS'
  ## The plasma somalogic surrogate signature for the healthy index
  SIGNATURE.IN.PATH = 'Classification/proteomic_surrogates/healthy.index.plasma.surrogate.somaId.RDS'
  ## The Ferrucci significance table
  TABLE.IN.PATH = 'Reference/ferrucci/raw/acel12799-sup-0004-TableS3.txt'
  
  SUPPLEMENTAL.FIGURE.5a.OUT.PATH = 'Paper_1_Figures/Supplemental_Figure_5/S5a.pdf'
  SUPPLEMENTAL.FIGURE.5b.OUT.PATH = 'Paper_1_Figures/Supplemental_Figure_5/S5b.pdf'
}

# Source utilities
source('scripts/util/Enrichment/hyperGeo.R')
source('scripts/util/Signatures/get_signature_scores.R')
# Figure 5a -- Forest plots of signature scores for each study in the meta-analysis for somalogic grey module proteins

## Load meta-analysis result and comparison group pairs
results = readRDS(META.ANALYSIS.Z.SCORE.IN.PATH)
cgps = readRDS(CGPS.IN.PATH)

## Get a map between the study name and its corresponding disease
### Create an empty vector to hold these names
studiess = c()
### For each disease
for(disease in names(cgps)) {
  ### Get the name of the studies for that disease
  studies = names(cgps[[disease]])
  ### Create a map from the study name to its corresponding disease
  new_studies = rep(disease, length(studies))
  names(new_studies) = studies
  ### Add the map to the empty vector
  studiess = c(studiess, new_studies)
}

## Get the sample size for each study
### For each disease
sizess = lapply(names(cgps), function(disease) {
  ### For each study of that disease
  studies = names(cgps[[disease]])
  ### Get the number of samples in that study
  sizes = sapply(cgps[[disease]], function(study) {
    length(unlist(study))
  })
  ### Name the sizes vector
  names(sizes) = names(cgps[[disease]])
  return(sizes)
})
sizess = unlist(sizess)

## Get the effect sizes and standard errors associated with each study
effects = results$datasetEffectSizes
ses = results$datasetEffectSizeStandardErrors

## Get the overall effect size and standard error
meta_effects = results$pooledResults[, 'effectSize', drop = FALSE]
meta_ses = results$pooledResults[, 'effectSizeStandardError', drop = FALSE]

## Instantiate a function to create the data frame used for the forest plot of a single feature's signature
get_df = function(feature) {
  ## Create an initial data frame with feature names, effect sizes, and standard error
  df1 = data.frame(study = colnames(effects), effect = effects[feature,], se = 1.96 * ses[feature,])
  
  ## Add the diseases associated with each study, the feature name, and the study size associated with the study
  df1$disease = factor(studiess[df1$study], c('DM1', 'MS', 'RA', 'sarcoid','summary',''))
  df1$feature = feature
  df1$study.size = sizess[df1$study]
  
  ## We also create a second data frame that is essentially a blank row to separate the diamond from the dots
  df2 = data.frame(study = '', effect = 0, 
                   se = 0, feature = feature,
                   disease = '', study.size = 0)
  
  
  ## We create a third data frame that just contains the meta_analysis effect size for display via a triange
  df3 = data.frame(study = 'Summary', effect = meta_effects[feature, ], 
                   se = 1.96 * meta_ses[feature, ], feature = feature,
                   disease = 'summary', study.size = 50)
  
  
  ## We put the data frames together
  df = rbind(df1, df2)
  df = rbind(df, df3)
  
  ## We put all the studies together into a factor
  #df$study = factor(df$study, levels = rev(unique(df$study)))
  return(df)
}

## We choose the features we want to show in the plot, and get their corresponding data frames for plotting
features = rownames(meta_effects)[c(2,7,13,15)]#grep('^somalogic\\.grey\\.', rownames(meta_effects), value = T)
dfs = lapply(features, get_df)

## We combine these dataframes
df = Reduce(rbind,dfs)

df = df %>%
  arrange(effect) %>%
  mutate(feature = gsub("somalogic\\.grey\\.", "", feature)) %>%
  mutate(feature = gsub("tbnks\\.", "", feature)) %>%
  mutate(feature = gsub("beta", "b", feature)) %>%
  mutate(feature = replace(feature, feature == "somalogic.modules.purple", "PM2")) %>%
  mutate(feature = replace(feature, feature == "healthy.index", "Immune Health Metric")) %>%
  mutate(feature = factor(feature, levels = unique(feature)))

## We rename the features for easier viewing
#levels(df$feature) = levels(df$feature) %>%
#  gsub('^somalogic\\.grey\\.', '', .) %>%
#  gsub('\\.', ' ', .)

## We manually create the standard ggplot colors
hues = seq(15, 375, length = 6)
colors = hcl(h = hues, l = 65, c = 100)[1:5]

# organize studies by disease
df <- df %>% arrange(desc(disease))
df$study <- factor(df$study,c("Summary",setdiff(unique(df$study),"Summary")))

## And create the forest plot
p = ggplot(df, aes(x = effect, y = study, color = disease)) +
  geom_point(aes(size = study.size, shape = disease), show.legend = T) +
  scale_shape_manual(values = c(16, 16, 16, 16, 18, 16)) + # 16 is for a circle and 18 a triangle
  geom_errorbarh(aes(xmin=effect-se, xmax=effect+se), height=0, show.legend = F, size = 1) +
  scale_color_manual(values = c(colors,'transparent')) + # We want the dot at 0 in the empty row to be transparent (we make it 0 to avoid the warnings from using an NA)
  xlab('Effect Size') + ylab('Study') + 
  theme_bw() + geom_vline(xintercept = 0, linetype = 'dashed') + facet_wrap(~feature, nrow = 1) + # We have a dashed line at 0 to represent no effect
  theme(axis.ticks.y = element_blank())

ggsave(SUPPLEMENTAL.FIGURE.5a.OUT.PATH, p, device = 'pdf', height = 5, width = 9)

# Figure 5b -- Barplot of Jamboree result enrichments for surrogate signatures
## We get the results of doing the transcriptional surrogate signature gene set enrichment among features included in the meta-analysis
results = readRDS(META.ANALYSIS.ENRICHMENTS.IN.PATH) 

## We get the pvalues associated with each transcriptional surrogate signature
pvals = sapply(results, function(result) {result$p.value}) 

## And convert to a negative log10 pvalue
negative.log10.pvals = -log10(pvals)


## We create a data frame with the feature names and negative log10 pvalues
df = data.frame(feature = names(pvals), negative.log10.pvals = negative.log10.pvals)
## We arrange the features by pvalue
df = df %>%
  arrange(pvals) %>%
  mutate(feature = gsub("somalogic\\.grey\\.", "", feature)) %>%
  mutate(feature = gsub("tbnks\\.", "", feature)) %>%
  mutate(feature = gsub("beta", "b", feature)) %>%
  mutate(feature = replace(feature, feature == "somalogic.modules.purple", "PM2")) %>%
  mutate(feature = replace(feature, feature == "healthy.index", "Immune Health Metric")) %>%
  mutate(feature = factor(feature, levels = unique(feature)))


source("scripts/util/Plotting/tbnk_featurename_replace.R")
levels(df$feature) <- replace_tbnk_names(levels(df$feature))
#levels(df$feature)

df <- df %>% filter(feature != "microarray.classifier")

## We create the bar plots
p = ggplot(df, aes(x = feature, y = negative.log10.pvals)) + geom_bar(stat = 'identity') + 
  coord_flip() + ylab('-log10(p)') + xlab('Feature') +
  theme_bw() + 
  theme(
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15)
  )
ggsave(SUPPLEMENTAL.FIGURE.5b.OUT.PATH, p, device = 'pdf', height = 6, width = 6)




