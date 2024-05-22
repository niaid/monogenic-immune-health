get_score <- function(x) {
  x = t(scale(t(x)))
  return (colMeans(x, na.rm=T))
}

calculate_signature_score <- function(expr.matrix,positive.contributors,negative.contributors) {
  sig.neg <- rownames(expr.matrix) %in% negative.contributors
  cat("Number of negative contributors represented in the dataset:",sum(sig.neg),"out of",length(negative.contributors),"\n")
  sig.pos <- rownames(expr.matrix) %in% positive.contributors
  cat("Number of positive contributors represented in the dataset:",sum(sig.pos),"out of",length(positive.contributors),"\n")
  sig_scores <- get_score(rbind(-expr.matrix[sig.neg,],expr.matrix[sig.pos,]))
  return(sig_scores)
}

generate_pseudo_sample_data <- function(expr.matrix,num.samples=10,spread.factor=10) {
  expr.var <- apply(expr.matrix,1,var)
  sample.list <- list()
  for (i in 1:ncol(expr.matrix)) {
    sample.mean <- expr.matrix[,i]
    pseudo.samples <- t(sapply(sample.mean, rnorm, n=num.samples, sd=spread.factor*expr.var))
    rownames(pseudo.samples) <- rownames(expr.matrix)
    colnames(pseudo.samples) <- paste0("Group",i,".",1:num.samples)
    sample.list[[i]] <- pseudo.samples
  }
  return(do.call(cbind,sample.list))
}