# A Unified Metric of Human Immune Health

This repository contains the code accompanying the article:
Rachel Sparks, Nicholas Rachmaninoff, William W. Lau, Dylan C. Hirsch, Neha Bansal, ..., John S. Tsang: A Unified Metric of Human Immune Health (Nature Medicine, 2024)


### Input Data

As the diseases in the study are extremely rare, data are deposited in dbGAP to protect patient confidentiality. Subject-level data, including microarray gene expressions, Somalogic-based serum protein measurements, clinical information, and demographics, are made available through dbGAP request under accession number phs002732.v1.p1 (https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002732.v1.p1). Aggregated mean gene count matrices, protein concentrations, and circulating cell frequencies by condition, as well as disease signatures and module membership are available publicly on https://panmonogenic.yale.edu/. Users can submit their own proteomic and/or transcriptional data for calculation of the IHM surrogate scores via this webtool.

### Instructions
The workflow to create all figures can be run with Snakemake and Singularity for increased reproducibility. 

Instructions on using Snakemake can be found here - https://snakemake.readthedocs.io/en/stable/

The singularity container to recreate the pipeline can be downloaded from sylabs.io with the following command.

```
singularity pull library://nrachman/default/monogenic:0.4
```

The definition file has also been provided in Snakemake/singularity_definition.

The entire workflow to go from raw data to figures/tables is described in the file called Snakemake/Snakefile.

An example script of starting the Snakemake workflow using the univa grid engine is provided in sm_call. This can be changed to use another high performance computing environment by changing sm_call and cluster_config.json. Alternatively, can be run in a single interactive session by calling `snakemake --use-singularity`. You will need at least 30GB of ram for some of the jobs. You can see which jobs are more resource intensive in cluster_config.json.

Code to create figures can be found in Snakemake/scripts/Paper_Figures_2024_03_15. Refer to Snakefile to find all necessary preprocessing code and input data.

Additionally, some code must be run as separate Rmarkdown files provided in the Revisions/ directory.

A file showing which scripts and outputs correspond to the different figures can be found in nm.submission.code.figure.mapping.xlsx

### Terms of Use

By using this software, you agree this software is to be used for research purposes only. Any presentation of data analysis using the software will acknowledge the software according to the guidelines below.

Primary author(s): Rachel Sparks, Nicholas Rachmaninoff, William W. Lau, Dylan C. Hirsch, Neha Bansal

Organizational contact information: John Tsang (john.tsang AT nih.gov)

Date of release: May 21, 2024

Version: 1.0

License details: see LICENSE file

Description: scripts used to generate figures and tables in the above publication

Usage instructions: See "Input Data" and "Instructions" sections from this readme.

Disclaimer:

A review of this code has been conducted, no critical errors exist, and to the best of the authors knowledge, there are no problematic file paths, no local system configuration details, and no passwords or keys included in this code. This open source software comes as is with absolutely no warranty.


