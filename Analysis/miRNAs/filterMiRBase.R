setwd("/home/I864741/TCC/miRNA_predict/Analysis/miRNAs")

###Create an array with miRNA names and its sequences
mirbase <- read.delim("filtered_mature.fa",stringsAsFactors=FALSE, sep=" ", header=FALSE)
mirbase.names <- mirbase[,6]
names(mirbase.names) <- mirbase[,1]