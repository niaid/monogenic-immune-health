library(variancePartition)
library(Biobase)

setwd("/Volumes/SG_DATA/PROJECTS/Monogenic_Project/")
jive <- readRDS("Integration_output/jive/sample/jive.rds")
eset <- readRDS("Data/Microarray/data_analysis_ready/eset_batch_filtered.rds")
varpart <- readRDS("Classification/variance_partitions/microarray_features.rds")

#select stable feats

array.stable.names <- 
  varpart@row.names[which(varpart$Patient > varpart$Residuals)]
eset <- eset[featureNames(eset) %in% array.stable.names]
  
#select pats
eset <- eset[, match(colnames(jive$data$array), eset$visit_id)]

#create mats
z.scale.mat <- as.matrix(Matrix::Diagonal(x = jive$scale$z.scale$array))
frob.scale.val <- jive$scale$frob.scale.val$array

z.center.mat <- matrix(rep(jive$scale$z.center$array, each =  ncol(eset)), 
                       nrow = nrow(eset), ncol = ncol(eset), byrow = TRUE)

#perform matrix multiplication and what not
array.reconstrctd <- t(t(jive$data$array) %*% z.scale.mat * frob.scale.val + t(z.center.mat))

#assess differences
smoothScatter(exprs(eset), array.reconstrctd)
cor(as.vector(exprs(eset)), as.vector(array.reconstrctd))

diff <- exprs(eset) - array.reconstrctd
sum(diff > .00001)
