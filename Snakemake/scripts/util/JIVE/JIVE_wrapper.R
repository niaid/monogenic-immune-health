jive_wrapper <- function(data.list, 
                         z.score = TRUE, frob.scale = TRUE, save.scale.info = TRUE,
                         pdat, save.pdat = TRUE,
                         method = "perm",
                         id.col){
  #This function wraps around the r.jiv::jive function using my implementation of the frobenius scaling and saves different scaling information
  
  # arguments:
  #   data.list is a list of data matrices that will go into jive. Features in rows. Samples in columns.
  #   z.score denotes whether features should be standardized to mean 0 and variance of 1 prior to scaling
  #   frob.scale determines whether the data matrices should be scaled by their frobenius norm so that all data matrices have the same amount of total variation
  #   pdat is the phenotypic data that corresponds to the samples
  #   save.pdat refers to whether you would like to include the pdat as part of the output list
  #   method refers to the way in which the number of JIVE PC's is chosen. See jive function for more detal

  # Value
	#   returns a list much like the output of the r.jive::jive function; however with different scaling information and potentially pdat included


  require(r.jive)
  
  #make sure samples are in each data matrix match the pdat 
  for(data.type in data.list){
    stopifnot(exists("id.col")) # bc of lazy eval
    stopifnot(identical(colnames(data.type), pdat[[id.col]]))
  }
  
  
  #scale prior to jive--------------------------------------------
  scale.info <- list()
  #Do standard feature z-scoring so that features with higher magnitude don't have higher weight
  if(z.score){
    data.list <- lapply(data.list, function(x) t(scale(t(x))) )
    if(save.scale.info){
      scale.info$z.scale <- lapply(data.list, function(x) attributes(x)$`scaled:scale`)
      scale.info$z.center <- lapply(data.list, function(x) attributes(x)$`scaled:center`)
    }
  }
  
  #The scaling in the package causes very small values that are potentially numerically unstable,
  #so I implement my own that does the same thing, but is multiplied by a constant
  if(frob.scale){
    data.list <- numerically_stable_frobenius_scale(data.list)
    if(save.scale.info){
      scale.info$frob.scale.val <- 
        lapply(data.list, function(x) attributes(x)$frob.scale.val)
    }
  }
  
  #run JIVE -----------------------------
  jive.out <- jive(data.list, scale = FALSE, center = FALSE, method = method)
  print("jive done")
  
  #add names to the joint and individual
  names(jive.out$joint) <- names(data.list)
  names(jive.out$individual) <- names(data.list)
  
  #add rownames and colnames to everything after performing jive
  for(data.type in c("data", "joint", "individual")){
    for(i in seq_along(data.list)){
      dimnames(jive.out[[data.type]][[i]]) <- dimnames(data.list[[i]])
    }
  }
  
  #update scaling information
  jive.out$scale <- scale.info
 
  #add the phenotypic data ----------------------------------
  if(save.pdat){
    jive.out$pdat <- pdat
  }
  
  jive.out
}

#------------------------------------------------------------------------
#Scaling
#------------------------------------------------------------------------

numerically_stable_frobenius_scale <- function(data){
  # The scaling in the package causes very small values that are potentially numerically unstable
  # If you print the jive function from the r.jive library, you will see how scaling is performed.
  # I scale the same way as the package; however, I normalize the scaling value to smallest frobenius norm
  
  # arguments: 
  #   data: a list of expression matrices, features X samples
  
  # Value
  #   a list of expression matrices that has been scaled and 
  #   with new scale value attributes added
  
  n <- c()
  #d <- c()
  for (i in 1:length(data)) {
    n[i] <- nrow(data[[i]]) * ncol(data[[i]])
    #d[i] <- nrow(data[[i]])
  }
  
  scaleValues <- sapply(data, function(x){
    norm(x, type = "f") * sqrt(sum(n))
  })
  
  #This is what I added to their scaling functionality
  scaleValues <- scaleValues/ min(scaleValues)
  
  for(i in 1:length(data)){
    data[[i]] <- data[[i]]/scaleValues[i]
    attributes(data[[i]])$frob.scale.val <- scaleValues[i]
  }
  
  return(data)
}
