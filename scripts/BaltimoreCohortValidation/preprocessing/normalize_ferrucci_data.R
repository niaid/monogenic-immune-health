rm(list=ls())

path_infile = paste0("/Users/candiajm/ACTIVE/SOMALOGIC_DATA/Ferruci_Serum/DATA/")
path_outfile = paste0("/Users/candiajm/ACTIVE/SOMALOGIC_DATA/Ferruci_Serum/RESULTS/")

norm = c("Raw","Hyb","Hyb.MedNorm","Hyb.MedNorm.Cal","Hyb.Cal","Hyb.Cal.MedNorm")
n_norm = length(norm)

dilution = c("0.005","1","40")
n_dilution = length(dilution)

seq_id_remove = c("2795-23","3590-8","5071-3","5118-74","5073-30") # taken down by Somalogic as of 12/2016.

dataset_label = c("CHI-18-003")
dataset = c("CHI-18-003.20180227")
n_dataset = length(dataset)

RFU_measured = NULL
metadata_sample = NULL
for (i_dataset in 1:n_dataset) {
  # we read somamer metadata
  tmp = t(as.matrix(read.table(paste0(path_infile,
                                      dataset[i_dataset],".adat_Somamers.txt"),header=F,sep="\t")))
  
  colnames(tmp) = tmp[1,]
  tmp = tmp[-1,]
  seq_id = strsplit(tmp[,1],"_")
  seq_id = matrix(unlist(seq_id),ncol=2,byrow=T)[,1] # to remove the sub-sequence info (after the underscore)
  tmp[,1] = seq_id
  somamer_remove = tmp[,1]%in%seq_id_remove
  tmp = tmp[!somamer_remove,]
  if (i_dataset==1) {
    metadata_analyte=tmp
    n_analyte=nrow(metadata_analyte)
  } else {
    cat("Dataset ",i_dataset,": Analyte integrity check passed = ",sum(metadata_analyte[,1]==tmp[,1])==n_analyte,"\n")
  }
  
  tmp = as.matrix(read.table(paste0(path_infile,
                                    dataset[i_dataset],".adat_RFU.txt"),header=F,sep="\t"))
  tmp = tmp[,!somamer_remove]
  RFU_tmp = tmp
  
  tmp = as.matrix(read.table(paste0(path_infile,
                                    dataset[i_dataset],".adat_Samples.txt"),header=T,sep="\t"))
  o = order(tmp[,"PlateId"])
  tmp = tmp[o,]
  RFU_measured = rbind(RFU_measured,RFU_tmp[o,])
  
  # Code below was used to analyze John Tsang's data:
  # "/Users/candiajm/ACTIVE/SOMALOGIC_DATA/Tsang"
  
  # To separate individual plates (in the case of multiplate adat files):
  #plate_label = unique(tmp[,"PlateId"])
  #n_plate = length(plate_label)
  #plate_label_data = rep(NA,nrow(tmp))
  #for (i_plate in 1:n_plate) {
  #    plate_label_data[tmp[,"PlateId"]==plate_label[i_plate]] = paste0(dataset_label[i_dataset],"_",i_plate)
  #}
  #tmp[,"PlateId"] = plate_label_data
  #if (ncol(tmp)==23) {
  #    tmp = tmp[,-23] # to fix discrepancies in column data
  #}
  metadata_sample = rbind(metadata_sample,tmp)
}
# we remove plate "CHI-18-003_2" (doesn't have any samples for the Tsang/Monogenic study)
#sample_remove = metadata_sample[,"PlateId"]=="CHI-18-003_2"
#metadata_sample = metadata_sample[!sample_remove,]
#RFU_measured = RFU_measured[!sample_remove,]

# we check what's in each plate
plate_label = unique(metadata_sample[,"PlateId"])
n_plate = length(plate_label)
SampleType_freq = vector("list",n_plate)
SampleMatrix_freq = vector("list",n_plate)
for (i_plate in 1:n_plate) {
  select = metadata_sample[,"PlateId"]==plate_label[i_plate]
  SampleType_freq[[i_plate]] = table(metadata_sample[select,"SampleType"])
  SampleMatrix_freq[[i_plate]] = table(metadata_sample[select,"SampleMatrix"])
}

#index = which(metadata_sample[,"SampleType"]=="Calibrator") # to check calibrators (5 x plate)
calibrator = "160385"

# we relabel bridge (to make it compatible -when merging later on)
index = which(metadata_sample[,"SampleType"]=="QC_CHI") # to check bridge samples (4 x plate)
metadata_sample[index,"SampleId"] = "QC_CHI"

RFU = vector("list",n_norm) # fluorescence intensity matrix

# 1. raw data
RFU[[1]] = RFU_measured

analyte_select_HCE = metadata_analyte[,"Type"]=="Hybridization Control Elution"
metadata_analyte_HCE = metadata_analyte[analyte_select_HCE,]
n_analyte_HCE = nrow(metadata_analyte_HCE)
RFU_HCE = RFU[[1]][,analyte_select_HCE]

# 2. hybridization control normalization
n_sample = nrow(RFU[[1]])
RFU[[2]] = RFU[[1]]
SF_Hyb = rep(NA,n_sample)
HCE_ratio = matrix(rep(NA,n_sample*n_analyte_HCE),ncol=n_analyte_HCE)
for (i_plate in 1:n_plate) {
  sample_select = metadata_sample[,"PlateId"]==plate_label[i_plate]
  hyb_intra_ref = apply(RFU_HCE[sample_select,],2,median)
  HCE_ratio_this_dataset = HCE_ratio[sample_select,]
  SF_Hyb_this_dataset = SF_Hyb[sample_select]
  RFU_this_dataset = RFU[[2]][sample_select,]
  for (i_sample in 1:sum(sample_select)) {
    HCE_ratio_this_dataset[i_sample,] = hyb_intra_ref/RFU_this_dataset[i_sample,analyte_select_HCE]
    SF_Hyb_this_dataset[i_sample] = median(HCE_ratio_this_dataset[i_sample,])
    RFU_this_dataset[i_sample,] = SF_Hyb_this_dataset[i_sample]*RFU_this_dataset[i_sample,]
  }
  HCE_ratio[sample_select,] = HCE_ratio_this_dataset
  SF_Hyb[sample_select] = SF_Hyb_this_dataset
  RFU[[2]][sample_select,] = RFU_this_dataset
}

# HERE: we remove samples that don't belong to the target study (we kept them to this point for the HCE normalization)
sample_labels = as.matrix(read.table(paste0(path_infile,"Tsang_17_samples.txt"),header=F,sep="\t"))[,1]
sample_select = !metadata_sample[,"SampleId"]%in%sample_labels
metadata_sample = metadata_sample[sample_select,]

RFU[[1]] = RFU[[1]][sample_select,!analyte_select_HCE]
RFU[[2]] = RFU[[2]][sample_select,!analyte_select_HCE]
metadata_analyte = metadata_analyte[!analyte_select_HCE,]
n_analyte = nrow(metadata_analyte)

# we write metadata to output
write(t(rbind(colnames(metadata_analyte),metadata_analyte)),ncol=ncol(metadata_analyte),file=paste0(path_outfile,"Somamers.txt"),sep="\t")
write(t(rbind(colnames(metadata_sample),metadata_sample)),ncol=ncol(metadata_sample),file=paste0(path_outfile,"Samples.txt"),sep="\t")

# 3. we add median normalization
sample_type_id = metadata_sample[,"SampleType"]
sample_type_id[sample_type_id=="QC"] = metadata_sample[sample_type_id=="QC","SampleId"]
RFU[[3]] = RFU[[2]]
for (i_plate in 1:n_plate) {
  sample_type = unique(sample_type_id[metadata_sample[,"PlateId"]==plate_label[i_plate]])
  n_sample_type = length(sample_type)
  for (i_sample_type in 1:n_sample_type) {
    select_sample = (metadata_sample[,"PlateId"]==plate_label[i_plate]) & (sample_type_id==sample_type[i_sample_type])
    for (i_dilution in 1:n_dilution) {
      select_somamer = metadata_analyte[,"Dilution"] == dilution[i_dilution]
      data = RFU[[3]][select_sample,select_somamer,drop=F]
      data2 = data
      for (i in 1:ncol(data2)) {
        data2[,i] = data2[,i]/median(data2[,i])
      }
      for (i in 1:nrow(data)) {
        data[i,] = data[i,]/median(data2[i,])
      }
      RFU[[3]][select_sample,select_somamer] = data
    }
  }
}

# 4. we add inter-plate calibration
RFU[[4]] = RFU[[3]]
reference_all_plates = apply(RFU[[4]][metadata_sample[,"SampleId"]==calibrator,],2,median)
for (i_plate in 1:n_plate) {
  sample_select = metadata_sample[,"PlateId"]==plate_label[i_plate]
  reference = apply(RFU[[4]][sample_select&(metadata_sample[,"SampleId"]==calibrator),],2,median)/reference_all_plates
  for (i_analyte in 1:n_analyte) {
    RFU[[4]][sample_select,i_analyte] = RFU[[4]][sample_select,i_analyte]/reference[i_analyte]
  }
}

# 5. we perform inter-plate calibration without median normalization (except for non-samples).
RFU[[5]] = RFU[[2]]
sample_select = metadata_sample[,"SampleType"]!="Sample"
RFU[[5]][sample_select,] = RFU[[3]][sample_select,]
reference_all_plates = apply(RFU[[5]][metadata_sample[,"SampleId"]==calibrator,],2,median)
for (i_plate in 1:n_plate) {
  sample_select = metadata_sample[,"PlateId"]==plate_label[i_plate]
  reference = apply(RFU[[5]][sample_select&(metadata_sample[,"SampleId"]==calibrator),],2,median)/reference_all_plates
  for (i_analyte in 1:n_analyte) {
    RFU[[5]][sample_select,i_analyte] = RFU[[5]][sample_select,i_analyte]/reference[i_analyte]
  }
}

# 6. we add multiplate median normalization on Samples.
RFU[[6]] = RFU[[5]]
select_sample = metadata_sample[,"SampleType"]=="Sample"
for (i_dilution in 1:n_dilution) {
  select_somamer = metadata_analyte[,"Dilution"] == dilution[i_dilution]
  data = RFU[[6]][select_sample,select_somamer,drop=F]
  data2 = data
  for (i in 1:ncol(data2)) {
    data2[,i] = data2[,i]/median(data2[,i])
  }
  for (i in 1:nrow(data)) {
    data[i,] = data[i,]/median(data2[i,])
  }
  RFU[[6]][select_sample,select_somamer] = data
}

# we write RFU data to output
for (i_norm in 1:n_norm) {
  write(t(RFU[[i_norm]]),ncol=ncol(RFU[[i_norm]]),
        file=paste0(path_outfile,norm[i_norm],"_RFU.txt"),sep="\t")
}
