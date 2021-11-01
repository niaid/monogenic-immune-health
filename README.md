# monogenic-immune-health

# Multiomics integration of 22 immune-mediated monogenic diseases reveals an emergent axis of immune health in humans

This repository contains the code accompanying the article:
Rachel Sparks, Dylan C. Hirsch, Nicholas Rachmaninoff, ..., John S. Tsang: Multiomics integration of 22 immune-mediated monogenic diseases reveals an emergent axis of immune health in humans (In review)


### Input Data

As the diseases in the study are extremely, data are being deposited in DBGAP to protect patient confidentiality. The DBGAP accession number will be provided once assigned.

### Instructions
The workflow to create all figures can be run with Snakemake and Singularity for increased reproducibility. 

Instructions on using Snakemake can be found here - https://snakemake.readthedocs.io/en/stable/

The singularity container to recreate the pipeline can be downloaded from sylabs.io with the following command.

```
singularity pull library://nrachman/default/monogenic:0.2
```

The entire workflow to go from raw data to figures/tables is described in the file called Snakefile.

An example script of starting the Snakemake workflow using the univa grid engine is provided in sm_call. This can be changed to use another high performance computing environment by changing sm_call and cluster_config.json. Alternatively, can be run in a single interactive session by calling `snakemake --use-singularity`

Code to create figures and supplementary tables can be found in scripts/Paper_Figures. Refer to Snakefile to find all necessary preprocessing code and input data.

### Terms of Use

By using this software, you agree this software is to be used for research purposes only. Any presentation of data analysis using the software will acknowledge the software according to the guidelines below.

Primary author(s): Rachel Sparks, Dylan C. Hirsch, Nicholas Rachmaninoff

Organizational contact information: John Tsang (john.tsang AT nih.gov)

Date of release: November 3, 2021

Version: 1.0

License details: see LICENSE file

Description: scripts used to generate figures and tables in the above publication

Usage instructions: See "Input Data" and "Instructions" sections from this readme.

Disclaimer:

A review of this code has been conducted, no critical errors exist, and to the best of the authors knowledge, there are no problematic file paths, no local system configuration details, and no passwords or keys included in this code. This open source software comes as is with absolutely no warranty.


