library(tidyverse)

if(!exists("snakemake")){
  setwd("../../..")
  source("scripts/util/paper/parse_snakemake.R")
  parse_snakemake(rule = "figure_4_meta_analysis_cgp_n_subj")
}

CGPS.CLEAN.IN.PATH <- snakemake@input[["clean_cgps"]]
CGPS.ORIG.IN.PATH <- snakemake@input[["orig_cgps"]]

TAB.OUT.PATH <- snakemake@output[[1]]

cgps_orig <- readRDS(CGPS.ORIG.IN.PATH)
cgps_clean <- readRDS(CGPS.CLEAN.IN.PATH)


dat <- lapply(cgps_clean, function(disease){
  lapply(disease, function(study){
    data.frame(case = length(study$case), control = length(study$control),
    case_samples = paste(study$case, collapse = " "), control_samples = paste(study$control, collapse = " "))
  }) %>% bind_rows(.id = "study_platform")
}) %>% bind_rows(.id = "disease")

cgp_orig_dat <- lapply(cgps_orig, function(disease){
  lapply(disease, function(cgp){
    as.data.frame(cgp$study.info)
  }) %>% bind_rows()
}) %>% bind_rows(.id = "disease")

cgp_orig_dat <- cgp_orig_dat %>% 
        mutate(study_platform = paste(study, platform, sep = "."))

rm_cgps <- c(
  "GSE9006-Diabetes_Mellitus,_Type_1-PBMC_newly diagnosed_paired with 1 month follow up::GSE9006-Healthy-PBMC_unpaired",
  "Jam_human_RA_GSE26554-JIA-PBMC::Jam_human_RA_GSE26554-Control-PBMC",
  "Jam_human_RA_GSE61281-Psoriatric_arthritis-Whole_blood::Cutaneous psoriasis without arthritis_GSE61281-Cutaneous_psoriasis_without_arthritis-Whole_blood",
  "Jam_Human_RA_JIA-PBMC::Jam_Human_RA_Controls-PBMC",
  "Jam_human_RA_GSE26554-Oligoarticular JIA-PBMC::Jam_human_RA_GSE26554-Control-PBMC",
  "Jam_Human_RA_JIA-PBMC::Jam_Human_RA_Controls-PBMC"
)

cgp_orig_dat <- cgp_orig_dat %>% filter(!name %in% rm_cgps)

cgp_orig_summ <- cgp_orig_dat %>%
        group_by(disease, study, platform, study_platform) %>%
        summarise(cgps = paste(name, collapse = "\t"))

dat <- dat %>%
        right_join(cgp_orig_summ)

note_dat <- tribble(
  ~study, ~note,
  "GSE21942", "GSM545843, GSM545845 were removed as these were technical replicates of other samples in the study",
  "GSE30210", "Removed additional replicates such that each individual only had one sample. Selected last sample chronologically",
  "GSE8650", "Removed additional replicates such that each individual only had one sample. Selected last sample chronologically. GSM214490 and GSM214492 were removed as they were believed to have unreliable diagnoses according to the original publication",
  "GSE15645", "Removed patients who were experiencing clinical remission of symptoms",
  "GSE42834", "Removed patients with non-active sarcoid"
)

dat <- dat %>% left_join(note_dat)

dat <- dat %>% filter(disease != "SLE")

write_csv(dat, TAB.OUT.PATH)


#FIG.OUT.PATH <- "Pipeline_out/Paper_1_Figures/Supplemental_Figure_4/n_subj_per_study_meta_analysis.pdf"
#dat_long <- dat %>% gather(key = sample_type, value = count, -c("disease", "study"))
#
#
#pdf(FIG.OUT.PATH)
#gridExtra::grid.table(dat, rows= rep("", nrow(dat)))
#dev.off()
#
#p <- ggplot(dat_long, aes(y = count, x = study)) +
#        geom_col(aes(fill = sample_type), position = "dodge") +
#        geom_text(aes(label = count)) +
#        ylab("n") +
#        coord_flip() +
#        facet_grid(disease ~ 1, scales = "free_y", space = "free_y") +
#        theme_bw() +
#        theme(strip.background.x = element_blank(), strip.text.x = element_blank())
#
#ggsave(plot = p, filename = FIG.OUT.PATH, width = 4, height = 5)
