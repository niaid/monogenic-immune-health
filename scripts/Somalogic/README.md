# Somalogic Overview:

## Somalogic Data Analysis Pipeline:
1. Preprocess the data (log transform, place in a summarized experiment object, and divide into training, testing, and qc)  `preprocessing/preprocess_somalogic.R`
2. Remove outlier samples from the Somalogic training data `preprocessing/filter_somalogic_training.R`
3. Form modules of protiens using WGCNA `module_creation/somalogic_wgcna.R`
