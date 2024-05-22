jaccard <- function(x, y){
  length(intersect(x, y)) / length(union(x, y))
}

jaccardMat <- function(geneset_list){
  out_mat <- matrix(nrow = length(geneset_list), ncol = length(geneset_list))
  rownames(out_mat) <- names(geneset_list)
  colnames(out_mat) <- names(geneset_list)

  diag(out_mat) <- 1
  
  for(i in seq_along(geneset_list)){
    for(j in seq_along(geneset_list)){
      if(i <= j){
        next()
      }
      out_mat[i, j] <- out_mat[j, i] <- jaccard(geneset_list[[i]], geneset_list[[j]])
    }
  }

  out_mat
}

jaccardMat2 <- function(geneset_list1, geneset_list2){
  out_mat <- matrix(nrow = length(geneset_list1), ncol = length(geneset_list2))
  rownames(out_mat) <- names(geneset_list1)
  colnames(out_mat) <- names(geneset_list2)

  for(i in seq_along(geneset_list1)){
    for(j in seq_along(geneset_list2)){
      out_mat[i, j]  <- jaccard(geneset_list1[[i]], geneset_list2[[j]])
    }
  }

  out_mat
}
