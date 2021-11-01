library(limma)
library(dplyr)

GENE.IN.PATH <- snakemake@input[["tm1"]]
PROTEIN.IN.PATH <- snakemake@input[["prot"]]

GENE.DAT.OUT.PATH <- snakemake@output[["tm1"]]
PROTEIN.DAT.OUT.PATH <- snakemake@output[["prot"]]

gene <- readRDS(GENE.IN.PATH)

protein <- readRDS(PROTEIN.IN.PATH)

cmat <- makeContrasts(
                      groupXCGD-groupSTAT1.GOF, 
                      group47CGD-groupSTAT1.GOF, 
                      (groupXCGD + group47CGD)/2 -groupSTAT1.GOF, 
                      levels = colnames(gene)

)

gene_cfit <- contrasts.fit(fit = gene, contrasts = cmat)
gene_cfit <- eBayes(gene_cfit)

protein_cfit = contrasts.fit(fit = protein, contrasts = cmat)
protein_cfit <- eBayes(protein_cfit)

topTable(gene_cfit, coef = 1)

xcgd_gene_dat <- topTable(gene_cfit, coef = 1, number = Inf)
forty7cgd_gene_dat <- topTable(gene_cfit, coef = 2, number = Inf)

xcgd_prot_dat <- topTable(protein_cfit, coef = 1, number = Inf)
forty7cgd_prot_dat <- topTable(protein_cfit, coef = 2, number = Inf)


#I.TAC and IFIT1, STAT1
#grep("TAC", prot_dat$Target, ignore.case = T, value = T)
#
#grep("ifi", prot_dat$Target, ignore.case = T, value = T)
#grep("ifn", prot_dat$Target, ignore.case = T, value = T)
#grep("p56", prot_dat$Target, ignore.case = T, value = T)
#grep("inter", prot_dat$Target, ignore.case = T, value = T)
#grep("isg", prot_dat$Target, ignore.case = T, value = T)
#grep("56", prot_dat$Target, ignore.case = T, value = T)
#
#grep("STAT1", prot_dat$Target, ignore.case = T, value = T)


gene_dat <- 
        list(`X-CGD - STAT1 GOF` = xcgd_gene_dat, 
             `p47-CGD - STAT1 GOF` = forty7cgd_gene_dat) %>%
        bind_rows(.id = "comparison") %>%
        filter(module_name == "red") %>%
        mutate(module_name = "TM1: Interferon")


prot_dat <- 
        list(`X-CGD - STAT1 GOF` = xcgd_prot_dat, 
             `p47-CGD - STAT1 GOF` = forty7cgd_prot_dat) %>%
        bind_rows(.id = "comparison") %>%
        filter(Target %in% c("I-TAC", "STAT1")) %>%
        select(-c(Units, Dilution, Type, Organism))

readr::write_csv(prot_dat, PROTEIN.DAT.OUT.PATH)
readr::write_csv(gene_dat, GENE.DAT.OUT.PATH)



