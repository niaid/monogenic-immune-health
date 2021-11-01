library(parallel)
library(igraph)

get_upper_tri <- function(mat){
  mat[upper.tri(mat)]
}

get_stat_mat <- function(mat, margin){
  
  min <- apply(mat, margin, min)
  max <- apply(mat, margin, max)
  median <- apply(mat, margin, median)
  mad <- apply(mat, margin, mad)
  mean <- apply(mat, margin, mean)
  sd <- apply(mat, margin, sd)
  cv <- sd / mean
  
  stat.mat <- cbind(min, max, median, mad, mean, sd, cv)
  
  return(stat.mat)
}

pairwise_dist_boot_summary <- function(ref.dist.mat = NULL, boot.dist.list){
  #have to remove stuff as I go because the memory becomes overloaded
  
  if(!is.null(ref.dist.mat)){
    #check now because I remove for memory efficiency
    dimnames_check(ref.dist.mat, boot.dist.list)
    }
  
  dist.vecs <- sapply(boot.dist.list, get_upper_tri)
  rm(boot.dist.list)
  
  out.mat <- get_stat_mat(dist.vecs, 1)
  
  if(!is.null(ref.dist.mat)){
    ref.dist.vec <- get_upper_tri(ref.dist.mat)
    out.mat <- cbind(out.mat, ref.dist.vec)
  }
  return(out.mat)
}

single_mse <- function(ref.vec, test.vec){
  n <- length(ref.vec)
  
  mean((test.vec - ref.vec) ^ 2) / n
}

pairwise_dist_mse <- function(ref.dist.mat, boot.dist.list){
  boot.dist.vecs <- lapply(boot.dist.list, get_upper_tri)
  
  ref.dist.vec <- get_upper_tri(ref.dist.mat)
  
  sapply(boot.dist.vecs, function(x){
    single_mse(ref.dist.vec, x)
  })
}

dimnames_check <- function(ref.dist, boot.dist.list){
  sapply(boot.dist.list, function(dist.mat){
    stopifnot( identical(rownames(dist.mat), rownames(ref.dist)) )
    stopifnot( identical(colnames(dist.mat), colnames(ref.dist)) )
  })
  return(0)
}

plot_stat_mat <- function(stat.mat, title){
  boot.cv <- stat.mat[, "cv"]
  boot.median <- stat.mat[,"median"] 
  ref.dist <- stat.mat[, "ref.dist.vec"]
  
  hist(boot.cv, main = title)
  smoothScatter(boot.median, ref.dist, main = title)
  smoothScatter(boot.median - ref.dist, ref.dist, main = title)
}

#---------------------------------------------------------------
#neighbor related
#---------------------------------------------------------------

neighbor_maintained <- function(ref.mst, mst.boot.list){
  
  #distance is set to NA so that every edge is distance 1
  ref.mst.dist <- distances(ref.mst, weights = NA)
  
  mst.boot.dists <- mclapply(mst.boot.list, mc.cores = detectCores(), function(graph.){
    dist.mat <- distances(graph., weights = NA)
    
    #make sure dimnames same
    stopifnot(identical(rownames(dist.mat), rownames(ref.mst.dist)))
    stopifnot(identical(colnames(dist.mat), colnames(ref.mst.dist)))
    return(dist.mat)
  })
  
  #select upper triangular
  ref.mst.dist <- get_upper_tri(ref.mst.dist)
  mst.boot.dists <- lapply(mst.boot.dists, get_upper_tri)
  
  #turn boot dist into matrix
  mst.boot.dists <- do.call(cbind, mst.boot.dists)
  
  #select neighbors in original graph
  ref.neighbor.indices <- which(ref.mst.dist == 1)
  
  #select ref neighbors in bootstrap iters
  mst.boot.dists <- mst.boot.dists[ref.neighbor.indices,]
  
  return(mst.boot.dists)
}



#---------------------------------------------------------------
#Not used, overly complicated ---------
#some bits of code may be useful for future functions though
#---------------------------------------------------------------

avg_ego_dist_single1 <- function(graph., ref.ego.list, nodes.ref){
  #does not currently work
  lapply(1:length(ref.ego.list), function(i){
    
    ego.member.indices <- get_vertex_indices(graph., ref.ego.list[[i]])
    node.index <- get_vertex_indices(graph., nodes.ref[i])
    return(list(node.index, ego.member.indices))
    #return(graph.)
    
    dists <- distances(graph., node.index, ego.member.indices, weights = NA)
    dists <- get_upper_tri(dists)
    avg.dist <- mean(dists)
    avg.dist
  })
}


get_vertex_indices <- function(g, vertex.names){
  #g is graph
  sapply(vertex.names, function(vertex.name){
    which(names(V(g)) == vertex.name)
  })
  
}


avg_ego_dist_single <- function(graph., ref.ego.list, nodes.ref){
  sapply(1:length(ref.ego.list), function(i){
    
    #Converting it to numeric gives the index I believe
    ego.member.indices <- as.numeric(ref.ego.list[[i]])
    node.index <- i
    #return(node.index)
    #return(ego.member.indices)
    #return(graph.)
    
    dists <- distances(graph., node.index, ego.member.indices, weights = NA)
    dists <- get_upper_tri(dists)
    avg.dist <- mean(dists)
    avg.dist
  })
}


check_vertex_names_same <- function(ref.mst, mst.boot.list){
  ref.vertex.names <- get.vertex.attribute(ref.mst, "name")
  
  boot.mst.vertex.names <- sapply(mst.boot.list, function(x){
    get.vertex.attribute(x, "name")
  })
  
  sum.not.equal.by.vertex <- sapply(1:length(ref.vertex.names), function(i){
    sum(boot.mst.vertex.names[i,] != ref.vertex.names[i])
  })
  
  stopifnot(sum(sum.not.equal.by.vertex) == 0)
}


ego_maintained <- function(ref.mst, mst.boot.list, ord = 1){
  #ord is order in the ego call- # of links away that is considered part of ego(i.e. neighborhood)
  
  #make sure vertex names are same across all
  check_vertex_names_same(ref.mst, mst.boot.list)
  
  #This would be used with avg_ego_dist_single1
  #nodes.ref <- V(ref.mst)
  #ref.ego.list <- ego(ref.mst, nodes = nodes.ref, order = ord)
  
  ref.ego.list <- ego(ref.mst, order = ord)
  
  mst.dists <- mclapply(mst.boot.list, mc.cores = detectCores(), function(g){
    avg_ego_dist_single(g, ref.ego.list, nodes.ref)
  })
  
  return(mst.dists)
  
}
