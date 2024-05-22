# Load libraries
library(ggplot2)
library(Biobase)
library(ggpubr)
library(dplyr)
library(tidyr)


# Load paths
if(exists("snakemake")){
  AGING.ESET.IN.PATH = snakemake@input[[1]]#'Reference/ferrucci/processed/aging_eset.RDS'
  PROTEOMIC.SIGNATURE.IN.PATH = snakemake@input[[2]]#'Classification/proteomic_surrogates/healthy.index.plasma.surrogate.somaId.RDS'
  META.ANALYSIS.Z.SCORE.IN.PATH = snakemake@input[[3]]#'Reference/jamboree/analysis_output/results/jamboree_z_score_results.RDS'
  CGPS.IN.PATH = snakemake@input[[4]]#'Reference/jamboree/data_analysis_ready/cgps_clean.RDS'
  
  FIGURE.5a.OUT.PATH = snakemake@output[[1]]#'Paper_1_Figures/Figure_5/5.a.pdf'
  FIGURE.5c.OUT.PATH = snakemake@output[[2]]#'Paper_1_Figures/Figure_5/5.c.pdf'
  FIGURE.5d.OUT.PATH = snakemake@output[[3]]#'Paper_1_Figures/Figure_5/5.d.pdf'
}

source('scripts/util/Signatures/get_signature_scores.R')
# Figure 5a -- Baltimore Healthy Aging Study, age versus HI surrogate signature scatterplot with regression line

## Load the Baltimore Healthy Aging Study eset and the somalogic Healthy Index plasma surrogate signature
eset = readRDS(AGING.ESET.IN.PATH)
healthy.index.surrogate.signature = readRDS(PROTEOMIC.SIGNATURE.IN.PATH)

## Subset the eset to just the samples (not QC, Calibrators, or buffers)
eset = eset[,eset$SampleType == 'Sample']

## Extract the scores from the signature
X = t(exprs(eset))
scores = util.get_signature_score(X, healthy.index.surrogate.signature)
ages = eset$Age

## Put the healthy index surrogate scores into a data frame
df = data.frame(Age = ages, Healthy.Index = scores)

## We create the plot
p = ggplot(df, aes(x = Age, y = Healthy.Index)) + geom_point() + 
  geom_smooth(method = 'lm', formula = y ~ x) + theme_bw() +
  stat_cor(label.x = 70, label.y = .55, label.sep = '\n', output.type = 'text') + 
  ylab('Proteomic Healthy Index Surrgoate') + xlab('Age')

## And save it
ggsave(FIGURE.5a.OUT.PATH, p, height = 4, width = 4)

# Figure 5b -- cartoon description of the meta-analysis

# Figure 5c -- Predictive ability on external data sets of transcriptional surrogate signatures for top predictors in healthy index
# (Note that the text box saying * p < .10, ** p < .05, *** p < .01 was made manually as it looks nicer with a proper text editor
# than via using ggplot)

## Load metaintegrator meta-analysis scores
results = readRDS(META.ANALYSIS.Z.SCORE.IN.PATH)
results = results$pooledResults # Get the list element holding the meta-analysis results

## Create a data frame with the feature names, effect sizes, a 95% confidence interval, and a BH-adjusted pvalue
df = data.frame(feature = rownames(results), 
                effect.size = results$effectSize, 
                se = 1.96 * results$effectSizeStandardError, 
                pval = results$effectSizeFDR)

## Remove the microarray classifier (which served as our "positive control" classifier to make sure that we weren't losing
## signal by creating surrogate signatures of other data type)
df = df %>% filter(feature != 'microarray.classifier')

## Manually clean up the feature names to be more clear and presentable
df = df %>% 
  mutate(feature = feature %>%
           gsub('somalogic\\.grey\\.','', .) %>%
           gsub('tbnks\\.','', .) %>%
           gsub('\\.',' ', .) %>%
           gsub('nk_cells_percent','NK Cells (%)', .) %>%
           gsub('nk_cells_abs','NK Cells (#)', .) %>%
           gsub('healthy index','Immune Health Metric', .) %>%
           gsub('somalogic modules purple','PM2', .) %>%
           gsub('beta','b', .) %>%
           gsub('PC1','jPC1', .) %>%
           gsub('rdw','RDW', .)) %>%
  mutate(feat_group = ifelse(feature %in% c("jPC1", "Immune Health Metric"), 0, 1))
    

## Order the features by effect size
df = df %>% 
  arrange(desc(effect.size)) %>%
  mutate(feature = factor(feature, levels = feature))

## Create columns for the number of stars to put next to each pvalue
df = df %>%
  mutate(p.value.stars = '') %>%
  mutate(p.value.stars = p.value.stars %>%
           replace(pval < .10,'*') %>%
           replace(pval < .05,'**') %>%
           replace(pval < .01,'***'))

levels(df$feature) <- levels(df$feature) %>%
        gsub(pattern = "", replacement = "")

## Create the plot
p = ggplot(df, aes(x = effect.size, y = feature, text = p.value.stars)) +
  geom_point() +
  geom_errorbarh(aes(xmin=effect.size-se, xmax=effect.size+se), height=0) +
  geom_text(aes(label = p.value.stars), nudge_y = .2, size = 5, show.legend = TRUE) + 
  ylab('Parameter') + xlab('Effect Size') +
  theme_bw() + scale_shape_manual(values = c(15,16,17,18)) +
  facet_grid(feat_group~1, scales = "free_y", space = "free") +
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_blank()
        )

## Save the  plot
ggsave(FIGURE.5c.OUT.PATH, p, height = 4, width = 6)
  
# Figure 5d -- Forest plots of signature scores for each study in the meta-analysis

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
  df$study = factor(df$study, levels = rev(levels(df$study)))
  return(df)
}

## We choose the features we want to show in the plot, and get their corresponding data frames for plotting
features = c("tbnks.nk_cells_abs","tbnks.rdw","somalogic.modules.purple","healthy.index")
dfs = lapply(features, get_df)

## We combine these dataframes
df = Reduce(rbind,dfs)
df$feature = factor(df$feature, features)

## We rename the features for easier viewing
levels(df$feature) = c('NK Cells (#)','RDW','PM2','Immune Health\nMetric')

## We manually create the standard ggplot colors
hues = seq(15, 375, length = 6)
colors = hcl(h = hues, l = 65, c = 100)[1:5]

## And create the forest plot
p = ggplot(df, aes(x = effect, y = study, color = disease)) +
  geom_point(aes(size = study.size, shape = disease), show.legend = T) +
  scale_shape_manual(values = c(16, 16, 16, 16, 18, 16)) + # 16 is for a circle and 18 a triangle
  geom_errorbarh(aes(xmin=effect-se, xmax=effect+se), height=0, show.legend = F, size = 1) +
  scale_color_manual(values = c(colors,'transparent')) + # We want the dot at 0 in the empty row to be transparent (we make it 0 to avoid the warnings from using an NA)
  xlab('Effect Size') + ylab('Study') + 
  theme_bw() + geom_vline(xintercept = 0, linetype = 'dashed') + facet_wrap(~feature, nrow = 1) + # We have a dashed line at 0 to represent no effect
  theme(axis.ticks.y = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))

## Save the  plot
ggsave(FIGURE.5d.OUT.PATH, p, height = 4, width = 9)
