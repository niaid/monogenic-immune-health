prcomp_list_varfilter <- function(prcomp.list, cutoff = .03){
  # This function takes a list of prcomp objects and filters the x matrix such that only PC's explaining 
  # a ceratin percentage of the variation are kept.
  
  #Arguments:
  # prcomp.list: a list of prcomp objects output
  # cutoff: the amount of variation that must be explained by a PC 
  #    for it to be kept in the x matrix
  
  #Value:
  # a list of prcomp objects where the x matrices have been filtered as described above
  
  
  lapply(prcomp.list, function(prcomp.obj) {
    #select the pc's that explain greater than 5% of the variation
    keep.pcs <- prcomp.obj$sdev^2 / sum(prcomp.obj$sdev^2) > cutoff
    prcomp.obj$x <- prcomp.obj$x[, keep.pcs, drop = FALSE]
    
    prcomp.obj
  })
}

