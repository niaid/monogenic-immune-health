## Extract sample data from files obtained
## internally from OMiCC

# Set globals
## CGP data
DATA.IN.PATHS = list(
  DM1 = snakemake@input[[1]],#'Reference/jamboree/raw/jam_human_DM1.txt',
  MS = snakemake@input[[2]],#'Reference/jamboree/raw/jam_human_MS.txt',
  RA = snakemake@input[[3]],#'Reference/jamboree/raw/jam_human_RA.txt',
  sarcoid = snakemake@input[[4]],#'Reference/jamboree/raw/jam_human_sarcoid.txt',
  SLE = snakemake@input[[5]]#'Reference/jamboree/raw/jam_human_SLE.txt'
)

## Compiled jamboree data
DATA.OUT.PATH = snakemake@output[[1]]#'Reference/jamboree/processed/jamboree_data.RDS'

## Instantiate a function to extract data from each file path and convert to a matrix
get_data = function(file.path) {
  data = read.csv(file.path, header = T, sep = "\t", comment.char = "!", row.names = 1)
  data = as.matrix(data)
}

## Apply this function over file paths
data = lapply(DATA.IN.PATHS, get_data)

## Save results
saveRDS(data, DATA.OUT.PATH)
