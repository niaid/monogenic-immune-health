
sig <- readRDS("/hpcdata/sg/sg_data/PROJECTS/Monogenic_Project/Classification/transcriptional_surrogates/surrogate_signatures.RDS")

HI <- sig$healthy.index

saveRDS(HI, "/hpcdata/sg/sg_data/PROJECTS/Monogenic_Project/Misc/healhy_index_transcriptomic_signature_2021_04_19.rds")
