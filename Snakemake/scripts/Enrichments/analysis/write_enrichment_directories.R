# Saves module enrichment results to tables

# Clear workspace
rm(list = ls())

# Make tables for the somalogic module gene set enrichments
base.path = 'Enrichments/somalogic_modules/'
es = readRDS('Data/Somalogic/analysis_output/enrichments/somalogic_module_gene_set_enrichments.RDS')
for(e in names(es)) {
  en = es[[e]]
  en$Set.Name = rownames(en)
  n = ncol(en)
  en = en[,c(n, 1:(n-1))]
  out.path = paste0(base.path, e, '.txt')
  write.table(en, row.names = FALSE, out.path, sep = '\t', quote = TRUE)
}

# Make tables for the somalogic module tissue set enrichments
base.path = 'Enrichments/somalogic_modules_tissue/'
es = readRDS('Data/Somalogic/analysis_output/enrichments/somalogic_module_tissue_set_enrichments.RDS')
for(e in names(es)) {
  en = es[[e]]
  en$Set.Name = rownames(en)
  n = ncol(en)
  en = en[,c(n, 1:(n-1))]
  out.path = paste0(base.path, e, '.txt')
  write.table(en, row.names = FALSE, out.path, sep = '\t', quote = TRUE)
}

# Make tables for the microarray module gene set enrichments
base.path = 'Enrichments/microarray_modules/'
es = readRDS('Data/Microarray/analysis_output/enrichments/microarray_module_gene_set_enrichments.RDS')
for(e in names(es)) {
  en = es[[e]]
  en$Set.Name = rownames(en)
  n = ncol(en)
  en = en[,c(n, 1:(n-1))]
  out.path = paste0(base.path, e, '.txt')
  write.table(en, row.names = FALSE, out.path, sep = '\t', quote = TRUE)
}
