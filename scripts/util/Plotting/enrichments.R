## Utility functions for plotting the top enrichments for a signature in a bar plot
require(ggplot2)

get_top_enrichments = function(enrichment, 
                              n = 3, 
                              gene.sets = c('btms', 'go', 'kegg', 'reactome')) {
  
  ## Gets the top enrichments from a list of enrichments for each gene set
  ##
  ## Inputs:
  ## enrichment - a data frame in the format output by the multiHyperGeoTest utility function in scripts/util/hyperGeo.R
  ## n - the maximum number of enrichments to pull for each signature from each gene set, and for each direction
  ## gene.sets - a character vector with the names of the gene set collections to pull from. They must match the names in the 'source' column
  ## in each data frame
  ##
  ## Outputs:
  ## A data.frame with the same columns as the original input, orgered by adjusted pvalue,
  ## but limited to features with an across-gene-set adjusted Pvalue of < .05.
  ## There are also columns 'term' with the names of the feature, NLP with the negative log 10 adjusted
  ## pvalue, and sign corresponding to the sign of the correlation
  
  # For each collection of gene sets
  dfs = lapply(gene.sets, function(gene.set) {
    
    # Get the enrichments for that gene set
    enrichment.subset = enrichment[enrichment$source == gene.set, , drop = F]
    
    # Sort the gene sets by adjust pvalue (adjusted across all gene sets)
    enrichment.subset = enrichment.subset[order(enrichment.subset$across.Adjusted.Pvalue, decreasing = FALSE), , drop = F]
    
    # If the number of enrichments for a gene set exceeds the maximum number we wish to take from each gene set
    if(nrow(enrichment.subset) >= n) {
      
      # Only take the top n gene set enrichments
      enrichment.subset = enrichment.subset[1:n, , drop = F]
    }
    
    return(enrichment.subset)
    
  })
  
  # Combine the results from all gene set collections
  df = Reduce(rbind, dfs)
  
  # If the results are empty
  if(nrow(df) == 0) {
    # Return NULL
    return(NULL)
  }
  
  # Otherwise, sort all the results by adjusted pvalues
  df = df[order(df$across.Adjusted.Pvalue, decreasing = TRUE), , drop = F]
  
  # And filter to just the gene sets with an adjusted pvalue of less than .05
  df = df[df$across.Adjusted.Pvalue < .05, , drop = F]
  
  # If there are no enrichmnets after filtering
  if(nrow(df) == 0) {
    # Return NULL
    return(NULL)
  }
  
  # Order the terms based on the order of the enrichments
  df$term = factor(rownames(df), levels = rownames(df))
  
  # Calculate the negative log 10 pvalues
  df$nlp = -log10(df$across.Adjusted.Pvalue)
  
  return(df)
}

get_top_enrichments_bidirectional = function(enrichments, 
                                             n = 3, 
                                             gene.sets = c('btms', 'go', 'kegg', 'reactome')) {
  
  ## Gets the top enrichments from a list of enrichments for each gene set
  ##
  ## Inputs:
  ## enrichments - a list with two data frames ('positive' and 'negative') as elements. Each data frame is in the format
  ## output by the multiHyperGeoTest utility function in scripts/util/hyperGeo.R. The positive data frame corresponds to
  ## the signature of genes positively correleated to the original feature, and the negative data frame corrsponds to the output
  ## by the signature of genes negatively correlated to the original feature
  ## n - the maximum number of enrichments to pull for each signature from each gene set, and for each direction
  ## gene.sets - a character vector with the names of the gene set collections to pull from. They must match the names in the 'source' column
  ## in each data frame
  ##
  ## Outputs:
  ## A data.frame with the same columns as the original input, orgered by adjusted pvalue,
  ## but limited to features with an across-gene-set adjusted Pvalue of < .05.
  ## There are also columns 'term' with the names of the feature, NLP with the negative log 10 adjusted
  ## pvalue, and sign corresponding to the sign of the correlation
  
  # For both the positive and negative signatures
  dfs = lapply(c('positive', 'negative'), function(direction) {
    
    # Get the top enrichments
    enrichment = enrichments[[direction]]
    df = get_top_enrichments(enrichment, n = n, gene.sets = gene.sets)
    
    # Append the direction of the signature
    df$direction = ifelse(direction == 'positive', 1, -1)
    return(df)
  })
  
  # Combine the positive and negative results
  df = Reduce(rbind, dfs)
  return(df)
}


make_enrichment_bar_plot = function(enrichment, 
                                    n = 3, 
                                    gene.sets = c('btms', 'go', 'kegg', 'reactome')) {
  
  ## Convenience Negative log10 pvalue bar plot function to plot enrichment pvalues
  ##
  ## Inputs:
  ## enrichment - a data frame in the format output by the multiHyperGeoTest utility function in scripts/util/hyperGeo.R
  ## n - the maximum number of enrichments to pull for each signature from each gene set, and for each direction
  ## gene.sets - a character vector with the names of the gene sets collections to pull from. They must match the names in the 'source' column
  ## in each data frame
  ##
  ## Outputs:
  ## ggplot barplot
  
  # Get the top gene set enrichments from an enrichments data frame
  df = get_top_enrichments(enrichment, n, gene.sets)
  
  # Put the gene sets in order by negative log 10 pvalue
  df = df %>% 
    arrange(nlp) %>%
    mutate(term = factor(term, levels = term)) %>%
    mutate(source = factor(source))
  
  # Make the bar plots
  p = ggplot(df, aes(x = term, y = nlp, fill = source)) + geom_bar(stat = 'identity') + 
    coord_flip() + ylab('Negative log10 q-value') + xlab('Gene Set') + 
    geom_hline(aes(yintercept = -log10(.05), linetype = 'FDR = .05'), color = 'gray', size = 1.5) +
    scale_linetype_manual(values = 'dashed') +
    guides(fill = guide_legend('Gene Set\nSource')) + theme_bw()
  return(p)
}

make_enrichment_bar_plot_bidirectional = function(enrichments, 
                                    n = 3, 
                                    gene.sets = c('btms', 'go', 'kegg', 'reactome')) {
  
  ## Convenience Negative log10 pvalue bar plot function to plot enrichment pvalues
  ##
  ## Inputs:
  ## enrichments - a list with two data frames ('positive' and 'negative') as elements. Each data frame is in the format
  ## output by the multiHyperGeoTest utility function in scripts/util/hyperGeo.R. The positive data frame corresponds to
  ## the signature of genes positively correleated to the original feature, and the negative data frame corrsponds to the output
  ## by the signature of genes negatively correlated to the original feature
  ## n - the maximum number of enrichments to pull for each signature from each gene set, and for each direction
  ## gene.sets - a character vector with the names of the gene sets collections to pull from. They must match the names in the 'source' column
  ## in each data frame
  ##
  ## Outputs:
  ## ggplot barplot
  
  # Get the top gene set enrichments from enrichments data frames made from positive and negative correlate signatures
  df = get_top_enrichments_bidirectional(enrichments, n, gene.sets)
  
  # Place a sign on the gene set enrichments
  df$nlp = df$direction * df$nlp
  
  # Put the gene sets in order by negative log 10 pvalue
  df = df %>% 
    arrange(sign(nlp), nlp) %>%
    mutate(term = factor(term, levels = term)) %>%
    mutate(source = factor(source))
  
  # Make the bar plots
  p = ggplot(df, aes(x = term, y = nlp, fill = source)) + geom_bar(stat = 'identity') + 
    coord_flip() + ylab('(Signed) Negative log10 q-value') + xlab('Gene Set') + 
    geom_hline(aes(yintercept = -log10(.05), linetype = 'FDR = .05'), color = 'gray', size = 1.5) +
    geom_hline(aes(yintercept = log10(.05), linetype = 'FDR = .05'), color = 'gray', size = 1.5) +
    scale_linetype_manual(values = 'dashed') +
    guides(fill = guide_legend('Gene Set\nSource')) + theme_bw()
  return(p)
}
