#!/usr/bin/bash

module load snakemake
#module load R/3.6.0-goolf-1.7.20 
source ~/.moduleload/R3.6.0_v2

snakemake -j 205 --cluster-config cluster_config.json --cluster "qsub -terse -l quick,avx2,mem_free={cluster.mem_free},h_vmem={cluster.h_vmem} {cluster.parallel_opts}" --rerun-incomplete -s Snakefile.py
