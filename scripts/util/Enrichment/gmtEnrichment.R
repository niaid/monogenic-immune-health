read.gmt = function(gmt.file) {
  
  con = file(gmt.file, open = "r")
  
  lines = readLines(con)
  lines = strsplit(lines, '\t')
  
  processes = sapply(lines, function(x) {x[1]})
  gene.sets = lapply(lines, function(x) {tail(x, -3)})
  
  rm(lines)
  close(con)
  
  names(gene.sets) = processes
  return(gene.sets)
  
}

make.set.list = function(named.vector) {
  
  set.names = unique(named.vector)
  set.list = lapply(set.names, function(set.name) {names(named.vector)[named.vector == set.name]})
  names(set.list) = set.names
  return(set.list)
  
}

reconcile.set.lists = function(sl1, sl2) {
  
  sl1.all = unlist(sl1)
  sl2.all = unlist(sl2)
  
  sl1 = lapply(sl1, function(s) {s[s %in% sl2.all]})
  sl2 = lapply(sl2, function(s) {s[s %in% sl1.all]})
  
  sl1 = sl1[sapply(sl1, length) > 0]
  sl2 = sl2[sapply(sl2, length) > 0]
  
  return(list(sl1, sl2))
  
}

do.fisher.tests = function(sl1, sl2){
  
  sls = reconcile.set.lists(sl1, sl2)
  sl1 = sls[[1]]
  sl2 = sls[[2]]
  
  background = unlist(sl1)

  results = sapply(sl1, function(s1) {sapply(sl2, function(s2) {do.fisher.test(s1, s2, background)})}, simplify = TRUE)
  return(results)
}

do.fisher.test = function(s1, s2, background) {
  
  in.s1 = background %in% s1
  in.s1 = as.numeric(in.s1)
  in.s1 = factor(in.s1, levels = c(0,1))
  
  in.s2 = background %in% s2
  in.s2 = as.numeric(in.s2)
  in.s2 = factor(in.s2, levels = c(0,1))
  
  contingency.tab = table(in.s1, in.s2) #rows first argument in table, columns are second argument
  result = fisher.test(contingency.tab, alternative = "greater")
  p.value = result$p.value
  return(p.value)
  
}

adjust.p.value.matrix = function(p.value.matrix, adjustment.method = p.adjust.methods) {
  
  p.vals = unlist(p.value.matrix) 
  p.adjs = p.adjust(p.vals, adjustment.method)
  p.adj.matrix = matrix(p.adjs, ncol = ncol(p.value.matrix), nrow = nrow(p.value.matrix))
  rownames(p.adj.matrix) = rownames(p.value.matrix)
  colnames(p.adj.matrix) = colnames(p.value.matrix)
  return(p.adj.matrix)
  
}

get.significant.terms = function(p.adj.matrix, alpha) {
  terms = rownames(p.adj.matrix)
  signifiant.terms = lapply(colnames(p.adj.matrix), function(col) {
    significant = p.adj.matrix[, col] < alpha
    significant.pvals = p.adj.matrix[significant, col]
    p.val.adj = significant.pvals
    })
  names(signifiant.terms) = colnames(p.adj.matrix)
  return(signifiant.terms)
}

run.somalogic.enrichment = function(modules.path, gmt.path, eset.path) {
  modules = readRDS(modules.path)
  modules = modules[modules != 'grey']
  eset = readRDS(eset.path)
  names(modules) = fData(eset)[names(modules), 'EntrezGeneSymbol']
  
  gene.sets.from.wgcna = make.set.list(modules) 
  gene.sets.from.database = read.gmt(gmt.path)
  
  p.val.matrix = do.fisher.tests(gene.sets.from.wgcna, gene.sets.from.database)
  p.adj.matrix = adjust.p.value.matrix(p.val.matrix)
  significant.terms = get.significant.terms(p.adj.matrix, .05)
  return(significant.terms)
}
