plot_pc_var_expl_scree <- function(prcomp.obj, n.pcs, title){
  var.explained <- prcomp.obj$sdev^2 / sum(prcomp.obj$sdev^2)
  var.explained <- var.explained[1:n.pcs]
  
  index <- 1:n.pcs
  
  plot(index, var.explained, main = title)
}