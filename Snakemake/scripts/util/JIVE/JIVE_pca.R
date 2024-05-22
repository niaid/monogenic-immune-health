get_jive_pca <- function(jive.obj){
  #Performs PCA on the various outputs from JIVE
  # returns a list of prcomp objects
  
  #Concatenate joint matrices and perform PCA. Keep only Number of PC's determined by rank of JIVE decomposition
  joint.concat <- t(do.call(rbind, jive.obj$joint))
  pca.joint <- prcomp(joint.concat, center = FALSE, scale. = FALSE, rank. = jive.obj$rankJ)
  pca.joint$sdev <- pca.joint$sdev[seq_len(jive.obj$rankJ)] #keep only sdev of non-zero var PC's
  
  n.data.types <- length(jive.obj$individual)
  pca.ind <- list()
  for(i in 1:n.data.types){
    A <- t(jive.obj$individual[[i]])
    PCA <- prcomp(A, center = FALSE, scale. = FALSE, rank. = jive.obj$rankA[[i]])
    PCA$sdev <- PCA$sdev[seq_len(jive.obj$rankA[[i]])] #keep only sdev of non-zero var PC's
    pca.ind[[i]] <- PCA
  }
  names(pca.ind) <- paste0(names(jive.obj$individual), ".ind")
  out.pca.list <- c(list(joint = pca.joint), pca.ind)
  
  return(out.pca.list)
}

