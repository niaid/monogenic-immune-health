Tissue Signatures
================
Dylan Hirsch
5/16/2019

R Markdown
----------

``` r
library(Biobase)
```

    ## Warning: package 'Biobase' was built under R version 3.5.1

    ## Loading required package: BiocGenerics

    ## Warning: package 'BiocGenerics' was built under R version 3.5.1

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind,
    ##     colMeans, colnames, colSums, dirname, do.call, duplicated,
    ##     eval, evalq, Filter, Find, get, grep, grepl, intersect,
    ##     is.unsorted, lapply, lengths, Map, mapply, match, mget, order,
    ##     paste, pmax, pmax.int, pmin, pmin.int, Position, rank, rbind,
    ##     Reduce, rowMeans, rownames, rowSums, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which, which.max,
    ##     which.min

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

``` r
rm(list = ls())
```

``` r
source('../../util/Enrichment/hyperGeo.R')
gene_sets = readRDS('../../../Gene_sets/processed/tissue_gene_sets.RDS')
modules = readRDS('../../../Data/Somalogic/analysis_output/wgcna_results/modules.rds')
```

Test 1: are tissue-specific markers enriched in the grey module compared to pan tissue signals?
-----------------------------------------------------------------------------------------------

``` r
hits = names(modules)[modules == 'grey']
universe = gene_sets$all$protein
gene_set = unique(unlist(gene_sets$general$protein))
hyperGeoTest(hits, universe, gene_set)
```

    ## $p.value
    ## [1] 1
    ## 
    ## $observed.hits
    ## [1] 0
    ## 
    ## $expected.hits
    ## [1] NaN
    ## 
    ## $odds.ratio
    ## [1] NaN

``` r
hits = names(modules)[modules == 'grey']
universe = gene_sets$all$protein
gene_set = unique(unlist(gene_sets$medium$protein))
hyperGeoTest(gene_set, universe, hits)
```

    ## $p.value
    ## [1] 1
    ## 
    ## $observed.hits
    ## [1] 0
    ## 
    ## $expected.hits
    ## [1] NaN
    ## 
    ## $odds.ratio
    ## [1] NaN

``` r
hits = names(modules)[modules == 'grey']
universe = gene_sets$all$protein
gene_set = unique(unlist(gene_sets$strict$protein))
hyperGeoTest(gene_set, universe, hits)
```

    ## $p.value
    ## [1] 1
    ## 
    ## $observed.hits
    ## [1] 0
    ## 
    ## $expected.hits
    ## [1] NaN
    ## 
    ## $odds.ratio
    ## [1] NaN
