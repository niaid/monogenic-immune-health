pre-filtering
================

This markdown file outputs expressionsets of the microarray data for the training, testing, and QC(CHI control) samples. Batch correction is performed using comBat; however, the data is also saved without batch correction in case the user would like the more raw data. Data is summarized to the gene level by selecting a single probeset per gene using the output of get\_anno\_pick\_probes.R. Genes are then filtered to remove supposed lowly expressed genes based on a visual inspection of the histogram of intensities. Lastly, genes are removed if they have high technical variance, i.e. they are more variable in the QC samples than they are in the patient samples. Lastly, the output expressionsets are saved for future use.

``` r
source("scripts_nick/util/Processing/averageTechnicalReplicates.R")
source("scripts_nick/util/Processing/averageRepeatSamples.R")
options(stringsAsFactors = FALSE)
```

Input paths
===========

``` r
#path to microarray expressionset
ESET.IN.PATH <- "Data/Microarray/raw/eset_rma_with_pData.rds"

#picked probes- those most correlated with first pc of all probes for gene
PICKED.PROBES.IN.PATH <- "Data/Microarray/probeset/output/picked_probes.txt"

#Feature data to be added to expressionset
FEATURE.DATA.IN.PATH <- "Data/Microarray/probeset/pre_downloaded_ann/probe_annotations_full.csv" 
```

### Output paths

``` r
#directory where filtered data to be saved
out.dir <- "Data/Microarray/data_analysis_ready"
if(!dir.exists(out.dir)) {dir.create(out.dir)}

TRAINING.SAMPLE.PATH <- "Data/Microarray/data_analysis_ready/eset_training_sample.rds"
TRAINING.SUBJECT.PATH <- "Data/Microarray/data_analysis_ready/eset_training_subject.rds"
QC.PATH <- "Data/Microarray/data_analysis_ready/eset_qc.rds"
VALIDATION.SAMPLE.PATH <- "Data/Microarray/data_analysis_ready/eset_validation_sample.rds" 
VALIDATION.SUBJECT.PATH <- "Data/Microarray/data_analysis_ready/eset_validation_subject.rds" 

BATCH.TRAINING.SAMPLE.PATH <- "Data/Microarray/data_analysis_ready/eset_batch_training_sample.rds"
BATCH.TRAINING.SUBJECT.PATH <- "Data/Microarray/data_analysis_ready/eset_batch_training_subject.rds"
BATCH.QC.PATH <- "Data/Microarray/data_analysis_ready/eset_batch_qc.rds"
BATCH.VALIDATION.SAMPLE.PATH <- "Data/Microarray/data_analysis_ready/eset_batch_validation_sample.rds"  
BATCH.VALIDATION.SUBJECT.PATH <- "Data/Microarray/data_analysis_ready/eset_batch_validation_subject.rds"  
```

``` r
eset <- readRDS(ESET.IN.PATH)
```

``` r
picked <- read.table(PICKED.PROBES.IN.PATH, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
```

``` r
fdata <- read.csv(FEATURE.DATA.IN.PATH, header = TRUE)
```

add feature data
================

``` r
fdata <- as(fdata[match(featureNames(eset), fdata$ID),], "AnnotatedDataFrame")
featureData(eset) <- fdata
```

Sample filtering
================

### Remove low rin-value samples

``` r
rin.cutoff <- 6
stripchart(eset$rin_value)
title("rin value")
abline(v = 6)
```

![](pre_filtering_files/figure-markdown_github/rin-filter-1.png)

``` r
keep.samples <- eset$rin_value > rin.cutoff | is.na(eset$rin_value)#some rin values are NA
eset <- eset[, keep.samples]
```

Detect outliers and flag for potential removal later
====================================================

Patient 54 appears odd, but comes from 2 separate visits so more likely biological

``` r
runPca <-function(expressionset, title, center = TRUE, scale. = TRUE){
  pca <- prcomp(as.matrix(t(exprs(expressionset))), 
                center = center, scale. = scale.)
  
  var.explained <- summary(pca)$importance[2,]
  qplot(pca$x[,1], pca$x[,2],
        geom = "text", 
        label = eset$patient_id,
        color = expressionset$condition) + 
    xlab(paste("PC1", var.explained[1])) +
    ylab(paste("PC2", var.explained[2])) + 
    ggtitle(title)
}

runPca(eset, title = "before cut")
```

![](pre_filtering_files/figure-markdown_github/pca-1.png)

### Samples above abline are flagged as outliers that can be removed in downstream analyses

``` r
#See WGCNA tutorial for motivation behind this
emat <- scale(t(exprs(eset)))
rownames(emat) <- make.names(eset$patient_id)
              
sampleTree = hclust(dist(emat), method = "average");
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2, cex = 0.4)

cuthi = 375
# Plot a line to show the cut
abline(h = cuthi, col = "red")
```

![](pre_filtering_files/figure-markdown_github/Sample-QC-hclust-outliers-1.png)

``` r
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = cuthi, minSize = 10)
table(clust)
```

    ## clust
    ##   0   1 
    ##   6 400

``` r
#Flag the outliers
hc.outlier <- (clust==0)
#add to pdata
pData(eset) <- pData(eset) %>%
  mutate(hc.outlier = hc.outlier)
```

Try out two batch correction schemes. Choose to use the combat correction
=========================================================================

``` r
batch <- factor(eset$assay_desc)
limma <- removeBatchEffect(exprs(eset), batch = batch)
combat <- ComBat(exprs(eset), batch = batch)
```

    ## Found14batches

    ## Adjusting for0covariate(s) or covariate level(s)

    ## Standardizing Data across genes

    ## Fitting L/S model and finding priors

    ## Finding parametric adjustments

    ## Adjusting the Data

``` r
diff.between.mats <- function(mat1, mat2, title = NULL){
  diff.mat <- mat1 - mat2
  
  diff.by.mean <- sapply(1:nrow(mat1), function(i){
    mean(abs(diff.mat[i,])) / mean(exprs(eset)[i,])
  })
  
  hist(diff.by.mean, main = title)
}
diff.between.mats(exprs(eset), limma, title = "limma vs original")
```

![](pre_filtering_files/figure-markdown_github/batch-correction-1.png)

``` r
diff.between.mats(exprs(eset), combat, title = "combat vs original")
```

![](pre_filtering_files/figure-markdown_github/batch-correction-2.png)

``` r
diff.between.mats(limma, combat, "limma vs combat")
```

![](pre_filtering_files/figure-markdown_github/batch-correction-3.png)

``` r
eset.batch <- new("ExpressionSet", exprs = combat)
phenoData(eset.batch) <- phenoData(eset)
featureData(eset.batch) <- featureData(eset)
```

Remove Patient 46 as we later learned that they had a medical procedure that could affect data
============================================================================

``` r
eset <- eset[, eset$patient_id != "P46"]
eset.batch <- eset.batch[, eset.batch$patient_id != "P46"]
```

Select probes most correlated with PC1 from pick probeset function
==================================================================

``` r
eset <- eset[rownames(exprs(eset)) %in% picked$ID, ]
eset.batch <- eset.batch[rownames(exprs(eset.batch)) %in% picked$ID, ]

picked <- picked[match(rownames(exprs(eset)), picked$ID),]

stopifnot(sum(rownames(exprs(eset)) != picked$ID) == 0)
stopifnot(sum(rownames(exprs(eset.batch)) != picked$ID) == 0)

featureNames(eset) <- picked$gene
featureNames(eset.batch) <- picked$gene
```

subset to training and make QC eset consisting of technical control samples
===========================================================================

``` r
eset.qc <- eset[, eset$analysis_group == "QC"]
eset.validation <- eset[, eset$analysis_group == "Validation"]
eset.train <- eset[, eset$analysis_group == "Discovery"]

eset.batch.qc <- eset.batch[, eset.batch$analysis_group == "QC"]
eset.batch.validation <- eset.batch[, eset.batch$analysis_group == "Validation"]
eset.batch.train <- eset.batch[, eset.batch$analysis_group == "Discovery"]

rm(eset)
rm(eset.batch)
```

Average the technical replicates
================================

``` r
# meta.cols1 <- setdiff(colnames(pData(eset.train)), c("visit_id","well"))
meta.cols1 = c('patient_id', 'assay_desc', 'condition', 'gender', "patient_age_at_time_of_blood_draw", 'race', 'visit_type')

eset.train <- averageTechnicalReplicates(eset.train, meta.cols = meta.cols1)
eset.validation <- 
  averageTechnicalReplicates(eset.validation, meta.cols = meta.cols1)

eset.batch.train <- 
  averageTechnicalReplicates(eset.batch.train, meta.cols = meta.cols1)
eset.batch.validation <- 
  averageTechnicalReplicates(eset.batch.validation, meta.cols = meta.cols1)
```

Change the sampleNames to the visit\_id
=======================================

``` r
sampleNames(eset.train) <- eset.train$visit_id
sampleNames(eset.validation) <- eset.validation$visit_id

sampleNames(eset.batch.train) <- eset.batch.train$visit_id
sampleNames(eset.batch.validation) <- eset.batch.validation$visit_id
```

Gene filtering
==============

### lowly expressed genes are genes to left of the redline of histograms

for rationale see:
<http://www.peterlangfelder.com/reducing-data-size-by-filtering-and-collapsing/> and <https://support.bioconductor.org/p/98635/> <https://f1000research.com/articles/5-1384/>

``` r
par(mfrow = c(2,2))
plotGeneQC <- function(expressionset, man_threshold, title){
  hist(apply(exprs(expressionset), 1, median), breaks = 250, main = title)
  abline(v = man_threshold, col = "coral4", lwd = 2, main = title)
  hist(exprs(expressionset), breaks = 250, main = title)
  abline(v = man_threshold, col = "coral4", lwd = 2, main = title)
}

#discovery samples
plotGeneQC(eset.train, man_threshold = 4.4, title = "Discovery Samples")
#Tempus control and QC CHI
plotGeneQC(eset.qc, man_threshold = 4.4, "Control Samples")
```

![](pre_filtering_files/figure-markdown_github/lowly-expressed-gene-filter-1.png)

``` r
#criterion to be filtered on later
lowly.expressed <- apply(exprs(eset.train), 1, median) < 4.4 # manual threshold determined from plots 
```

Remove genes with high technical variance
=========================================

### Genes to the left of line have high technical variance

We define high technical variance as those genes with higher variance in the control samples than in the across all discovery samples. This likely also removes genes showing significant batch effects as the controls were run across the different batches.

"Control" samples are a single healthy control that were from the Tsang lab. "CHI Control" are a pooled set of whole blood that are used for controls in microarray experiments at the Center for Human Immunology of the NIH.

To be considered highly technical variable, the genes must be more variable in BOTH of these controls

``` r
smoothScatter(apply(exprs(eset.train), 1, sd), 
              apply(exprs(eset.qc[,eset.qc$condition == "Control"]), 1, sd), 
              xlab = "sd in discovery samples", ylab = "sd in Tsang lab 'Control' samples")
abline(0, 1, col = "red")
```

![](pre_filtering_files/figure-markdown_github/plot-technical-var-1.png)

``` r
smoothScatter(apply(exprs(eset.train), 1, sd), 
              apply(exprs(eset.qc[,eset.qc$condition == "CHI Control"]), 1, sd), 
              xlab = "sd in discovery samples", ylab = "sd in 'CHI Control' samples")
abline(0, 1, col = "red")
```

![](pre_filtering_files/figure-markdown_github/plot-technical-var-2.png)

``` r
smoothScatter(apply(exprs(eset.qc[,eset.qc$condition == "Control"]), 1, sd), 
              apply(exprs(eset.qc[,eset.qc$condition == "CHI Control"]), 1, sd), 
              xlab = "sd in Tsang lab 'Control'", ylab = "sd in chi control")
abline(0,1)
```

![](pre_filtering_files/figure-markdown_github/plot-technical-var-3.png)

### filtering based to remove genes that are lowly expressed and have high technical variance

``` r
hi.technical.variance <- 
  apply(exprs(eset.train), 1, sd) < apply(exprs(eset.qc[,eset.qc$condition == "Control"]), 1, sd) & 
  apply(exprs(eset.train), 1, sd) < apply(exprs(eset.qc[,eset.qc$condition == "CHI Control"]), 1, sd) 

eset.train <- eset.train[!lowly.expressed & !hi.technical.variance, ]
eset.qc <- eset.qc[!lowly.expressed & !hi.technical.variance, ]
eset.validation <- eset.validation[!lowly.expressed & !hi.technical.variance, ]

eset.batch.train <- eset.batch.train[!lowly.expressed & !hi.technical.variance, ]
eset.batch.qc <- eset.batch.qc[!lowly.expressed & !hi.technical.variance, ]
eset.batch.validation <- eset.batch.validation[!lowly.expressed & !hi.technical.variance, ]
```

Average repeat samples in the training set
==========================================

``` r
#meta.cols2 <- setdiff(colnames(pData(eset.train)), 
#                     c("patient_id", "visit_id",
#                       "patient_age_at_time_of_blood_draw"))

meta.cols2 = c('patient_id', 'assay_desc', 'condition', 'gender', 'race')

eset.train.sample <- eset.train
eset.train.subject <- averageRepeatSamples(eset.train, 
                                           patient.id.col = "patient_id",
                                           meta.cols = meta.cols2)

eset.batch.train.sample <- eset.batch.train
eset.batch.train.subject <- averageRepeatSamples(eset.batch.train, 
                                                 patient.id.col = "patient_id",
                                                 meta.cols = meta.cols2)

eset.validation.sample <- eset.validation
eset.validation.subject <- averageRepeatSamples(eset.validation, 
                                                 patient.id.col = "patient_id",
                                                 meta.cols = meta.cols2)

eset.batch.validation.sample <- eset.batch.validation
eset.batch.validation.subject <- averageRepeatSamples(eset.batch.validation, 
                                                 patient.id.col = "patient_id",
                                                 meta.cols = meta.cols2)
```

### save these expressionsets

``` r
saveRDS(eset.train.sample, file = TRAINING.SAMPLE.PATH)
saveRDS(eset.train.subject, file = TRAINING.SUBJECT.PATH)
saveRDS(eset.qc, file = QC.PATH)
saveRDS(eset.validation.sample, file = VALIDATION.SAMPLE.PATH)
saveRDS(eset.validation.subject, file = VALIDATION.SUBJECT.PATH)

saveRDS(eset.batch.train.sample, file = BATCH.TRAINING.SAMPLE.PATH)
saveRDS(eset.batch.train.subject, file = BATCH.TRAINING.SUBJECT.PATH)
saveRDS(eset.batch.qc, file = BATCH.QC.PATH)
saveRDS(eset.batch.validation.sample, file = BATCH.VALIDATION.SAMPLE.PATH)
saveRDS(eset.batch.validation.subject, file = BATCH.VALIDATION.SUBJECT.PATH)
```
