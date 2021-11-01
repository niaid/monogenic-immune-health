#This is how I generated the summary
#snakemake -n -D > workflow_summary.txt

library(tidyverse)

wf <- read_tsv("../workflow_summary.txt", skip = 1)

wf %>% as.data.frame() %>% head()


