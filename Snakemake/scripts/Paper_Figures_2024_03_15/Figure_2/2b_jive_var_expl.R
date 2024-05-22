library(tidyverse)
library(r.jive)

JIVE.PC.PATH <- snakemake@input[["jive_pcs"]]#"Integration_output/jive/subject/prcomp_list.rds"
JIVE.PATH <- snakemake@input[["jive"]]#"Integration_output/jive/subject/jive.rds"

FIG.OUT.PATH <- snakemake@output[[1]]#"Paper_1_Figures/Figure_3/jive_var_exp.pdf"

prcomp.list <- readRDS(JIVE.PC.PATH)
jive <- readRDS(JIVE.PATH)

#std dev of each pc
prcomp.list$joint$sdev

#The variance explained plot as done by the package
#showVarExplained(jive)
#Calculating the variance explained
#see r.jive::showVarExplained
result <- jive
l <- length(result$data)
VarJoint = rep(0, l)
for (i in 1:l) VarJoint[i] = norm(result$joint[[i]], type = "F")^2/norm(result$data[[i]], 
    type = "F")^2
VarIndiv = rep(0, l)
for (i in 1:l) VarIndiv[i] = norm(result$individual[[i]],
    type = "F")^2/norm(result$data[[i]], type = "F")^2
VarResid = 1 - VarJoint - VarIndiv


#Put variance explained into data frame


dat <- data.frame(Joint = VarJoint, Individual = VarIndiv, Residual = VarResid,
                  data.type = c("WB Transcriptome", "Serum Proteins")) %>% 
  gather(key = component, value = var.explained, -data.type)

dat$label <- vector("character", nrow(dat))

#make label that contains both joint, individual residual and microarray/somalogic
for(i in seq_len(nrow(dat))){
  if(dat$component[[i]] != "Joint"){
    print(i)
    dat$label[[i]] <- paste(dat$data.type[[i]], dat$component[[i]])
  }else{
    dat$label[[i]] <- as.character(dat$component)[[i]]
  }
}

dat$label <- factor(dat$label, 
                    levels = c( "Serum Proteins Residual", "Serum Proteins Individual", "WB Transcriptome Residual", "WB Transcriptome Individual", "Joint"))

#The plot that goes into the figure
pdf(FIG.OUT.PATH, height =1.5, width = 6)
ggplot(dat, aes(x = data.type, y = var.explained, fill = label)) +
  geom_col() +
  theme_classic() +
  coord_flip() +
  scale_fill_manual(values = c("WB Transcriptome Residual" = "lightblue2", 
                               "WB Transcriptome Individual" = "skyblue3",
                               "Serum Proteins Residual" = "coral1",
                               "Serum Proteins Individual" = "coral3",
                               "Joint" = "lightsteelblue4")) +
ylab("Variance\nExplained")# +
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()
