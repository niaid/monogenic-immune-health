## Extract comparison group pair (CGP) info from files obtained
## internally from OMiCC

# Set globals
## Paths to OMiCC files
CGP.IN.PATHS = list(
  DM1 = snakemake@input[[1]],#'Reference/jamboree/raw/jam_human_DM1.txt',
  MS = snakemake@input[[2]],#'Reference/jamboree/raw/jam_human_MS.txt',
  RA = snakemake@input[[3]],#'Reference/jamboree/raw/jam_human_RA.txt',
  sarcoid = snakemake@input[[4]],#'Reference/jamboree/raw/jam_human_sarcoid.txt',
  SLE = snakemake@input[[5]]#'Reference/jamboree/raw/jam_human_SLE.txt'
)

## Compiled jamboree cgps
CGP.OUT.PATH = snakemake@output[[1]]#'Reference/jamboree/processed/jamboree_cgps.RDS'

# Helper function for getting CGP info
get_info = function(indices, keys, values) {
  ## Parse a set of keys value pairs, limited only to
  ## the indices given by 'indices'
  
  ## Get all keys of the proper indices
  keys = keys[indices]
  ## Get all values of the proper indices
  values = values[indices]
  ## Extract the key name: this is the last entry after a '_'
  keys = sapply(strsplit(keys,'_'), function(x) {x[[length(x)]]})
  
  ## Create a key-value pair map
  values = as.list(values)
  names(values) = make.names(keys)
  return(values)
}

# Instantiate a function to get info a associated with a CGP
extract_CGP = function(cgp.meta, cgp.number) {
  ## Collects all information on a certain cgp number from metadata
  ## Get the header corresponding to this CGP number
  header = paste0("^!CGP_", cgp.number, "_")
  
  ## Get all lines in the meta data corresponding to this CGP
  cgp.meta = grep(header, cgp.meta, value=T)
  ## Remove the CGP header from each line
  cgp.meta = gsub(header, "", cgp.meta)
  ## Split the rest of the line into key value pairs
  cgp.meta = strsplit(cgp.meta, '\t')
  
  ## Get all of the keys
  keys = sapply(cgp.meta, function(x) {x[[1]]})
  ## Get all of the associated values
  values = sapply(cgp.meta, function(x) {x[[2]]})
  
  ## Get all lines that are general metadata for the CGPs
  ## These are those in which the keys do not contain the keyword 'condition'
  index = !grepl('condition', keys)
  ## Parse these lines
  study.info = get_info(index, keys, values)
  
  ## Get all lines that are metadata for the cases
  ## These are those in which the keys contain the keyword 'condition1_sample_'
  index = grepl('condition1_sample_', keys)
  ## Parse these lines
  case.info = get_info(index, keys, values)
  
  ## Get all lines that are metadata for the controls
  ## These are those in which the keys contain the keyword 'condition1_sample_'
  index = grepl('condition2_sample_', keys)
  ## Parse these lines
  control.info = get_info(index, keys, values)
  
  ## Get all GSMs for the cases
  ## These correspond to the value in which the key is condition1_sample
  case.names = values[keys == 'condition1_sample']
  ## Get all GSMs for the controls
  ## These correspond to the value in which the key is condition2_sample
  control.names = values[keys == 'condition2_sample']
  
  ## Put all the metadata into a list
  cgp = list(study.info = study.info,
             case.info = case.info,
             control.info = control.info,
             case.names = case.names,
             control.names = control.names)
  
  return(cgp)
}

# Instantiate function to extract CGP info from a file path
extract_CGPs = function(file.path) {
  ## Search for all lines containing '!' as the first character;
  ## these are the lines with CGP meta data.
  cgp.meta = grep("^!", readLines(file.path), value=T)
  ## Get the number of CGPs in the file, given by the second entry in the first line of metadata
  num.cgp = as.numeric(unlist(strsplit(cgp.meta[1],"\t"))[2])
  ## For each CGP, extract all information on that cgp from the metadata. Note that the
  ## CGP numbers start at 0
  cgps = lapply(0:(num.cgp-1), function(cgp.number) {extract_CGP(cgp.meta, cgp.number)})
  ## Format the CGP names, and make them start at 1 rather than 0
  names(cgps) = paste0('CGP_', 1:num.cgp)
  return(cgps)
}

# Extract CGP info
cgps = lapply(CGP.IN.PATHS, extract_CGPs)

# Save results
saveRDS(cgps, CGP.OUT.PATH)
