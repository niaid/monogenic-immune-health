perm_2group <- function(x, y, n.perm, summary.stat){
  x.summ <- summary.stat(x)
  y.summ <- summary.stat(y)
  summ.diff <- x.summ - y.summ
  
  perm.diffs <- numeric(n.perm)
  xy <- c(x, y)
  x.len <- length(x)
  for(i in seq_along(perm.diffs)){
    xy.perm <- sample(xy, length(xy), replace = FALSE)
    
    x.perm <- xy.perm[seq_len(x.len)]
    x.perm.summ <- summary.stat(x.perm)
    
    y.perm <- xy.perm[setdiff(seq_len(length(xy.perm)), seq_len(x.len))]
    y.perm.summ <- summary.stat(y.perm)
    
    perm.diff.single <- x.perm.summ - y.perm.summ 
    perm.diffs[i] <- perm.diff.single
  }
  
  two.sided.p <- two_sided_p(perm.diffs, summ.diff)
  
  results <- list(summ.diff = summ.diff,
                  two.sided.p = two.sided.p)
}

two_sided_p <- function(perm.diffs, summ.diff){
  two.sided.num <- 2 * min((sum(perm.diffs <= summ.diff) + 1),
                           (sum(perm.diffs >= summ.diff) + 1))
  two.sided.denom <- length(perm.diffs) + 1
  two.sided.p <- two.sided.num / two.sided.denom
  two.sided.p
}