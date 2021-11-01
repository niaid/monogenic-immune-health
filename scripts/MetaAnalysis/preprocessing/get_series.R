## Parses and organizes series matrix files downloaded from GEO using the GEOquery library

# Set Globals
## The series matrix files from GEO used for the meta-analysis
GEO.IN.PATHS = snakemake@input[1:(length(snakemake@input)-1)] #Reference/jamboree/raw/series_matrices/
## The comparison groups pairs from the jamboree
CGPS.IN.PATH = snakemake@input[[length(snakemake@input)]]#'Reference/jamboree/processed/jamboree_cgps.RDS'

## The series matrices
SERIES.OUT.PATH = snakemake@output[[1]]#'Reference/jamboree/processed/series_matrix_list.rds'

# Load data
cgps <- readRDS(CGPS.IN.PATH)

# Initialize an empty list of series matrices
series.list <- list()

# For each disease
for(nm1 in names(cgps)) {
  # For each CGP under that disease
  for(nm2 in names(cgps[[nm1]])) {
    # Get the GEO study id associated with the CGP
    geo.id <- cgps[[nm1]][[nm2]]$study.info$study
    # Get the GEO platform associated with the CGP
    geo.platform <- cgps[[nm1]][[nm2]]$study.info$platform
    
    # Get the path for the series matrix from the snakemake inputs
    in.name <- paste(geo.id, geo.platform, 'series_matrix.txt', sep = '_')
    in.path <- grep(paste0('*\\/', in.name), GEO.IN.PATHS, value = T)
    
    # Read in the series matrix file
    series.matrix = read.table(in.path, header = TRUE, quote = '', comment.char = '', sep = '\t',
                               stringsAsFactors = FALSE)
    
    # Add the series matrix to the list
    series.name = paste(geo.id, geo.platform, sep = '.')
    series.list[[series.name]] <- series.matrix
  }
}

saveRDS(series.list, file = SERIES.OUT.PATH)

