setwd("D:/TCC/miRNA_predict/Analysis")

#RNAtable <- read.delim("D:/TCC/miRNA_predict/Analysis/RNAs/names-sequences.txt",stringsAsFactors = FALSE,sep=",")
#miRNAtable <- read.delim("D:/TCC/miRNA_predict/Analysis/miRNAs/filtered_mature.fa",stringsAsFactors = FALSE,sep=" ")
load("D:/TCC/miRNA_predict/Analysis/miRNAs/processInteractionData_output/preprocessed_data_mirTarBasev6.1_TarBasev7.RData")

write.table(tarbase.negative,file="tarbase-negative.csv",sep=' ')

write.table(tarbase.positive,file="tarbase-positive.csv",sep=' ')