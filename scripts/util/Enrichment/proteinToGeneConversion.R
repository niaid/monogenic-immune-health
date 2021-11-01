make_protein_to_gene_map = function(eset) {
  ## Makes a map between proteins and gene names, removing proteins that
  ## correspond to multiple genes.
  ## Inputs:
  ## eset - an ExpressionSet of somamer data, with featureData on each somamer
  ## Outputs:
  ## named character vector - a map from protein names to gene names
  
  # Extract the somamer information from the proteins eset
  feature.meta = fData(eset)
  
  # Get the entrez gene symbol for each protein
  protein.to.gene.map = feature.meta$EntrezGeneSymbol
  
  # Map these genes to the protein name
  names(protein.to.gene.map) = rownames(feature.meta)
  
  # Remove proteins corresponding to multiple genes
  genes = unique(protein.to.gene.map)
  genes.to.remove = union(grep(' ', genes, value = T), grep(',', genes, value = T))
  protein.to.gene.map = protein.to.gene.map[! protein.to.gene.map %in% genes.to.remove]
  
  return(protein.to.gene.map)
}

get_coherent_genes = function(protein.to.gene.map, module) {
  ## Get each gene for which all proteins that map to it are either
  ## in a protein module or out of it
  ## Inputs:
  ## protein.to.gene.map - named character vector as output from make_protein_to_gene_map
  ## module - character vector of protein names (as in the rownames of the eset used to derive
  ##          the protein to gene map) that belong to a protein module
  ## Outputs:
  ## character vector - the name of all genes (as specified in the fData of the eset used to derive
  ##                    the protein to gene map) for the proteins that map to the gene are either all
  ##                    in the module or all out of the module
  
  # Get the proteins in a module
  in.proteins = module[module %in% names(protein.to.gene.map)]
  
  # Get the proteins in the universe outside the module
  out.proteins = setdiff(names(protein.to.gene.map), in.proteins)
  
  # Get the genes that are mapped to by proteins in the module
  in.genes = protein.to.gene.map[in.proteins]
  
  # Get the genes that are mapped to by proteins outside of the module
  out.genes = protein.to.gene.map[out.proteins]
  
  # Get the symmetric differene to obtain genes for which either all corresponding proteins
  # are in the module or all corresponding proteins are not in the module
  coherent_genes = setdiff(protein.to.gene.map, intersect(in.genes, out.genes))
  
  return(coherent_genes)
}
