setwd("D:/TCC/miRNA_predict/Analysis/RNAs")

namesTable <- read.delim("D:/TCC/miRNA_predict/Analysis/RNAs/mart_names.txt",stringsAsFactors = FALSE,sep=",")
sequencesTable <- read.delim("D:/TCC/miRNA_predict/Analysis/RNAs/mart_sequences.txt",stringsAsFactors = FALSE,sep=",",header=FALSE, col.names=c("Gene stable ID","Sequence"))

mergedTable <- merge(x = namesTable, y = sequencesTable, by = "Gene.stable.ID")
cleanedTable <- mergedTable[(mergedTable$Sequence!='Sequence unavailable'),]
cleanedTable <- unique(cleanedTable)

save(cleanedTable,file="names-sequences.RData")
write.table(cleanedTable,file="names-sequences.txt",row.names=FALSE,sep=',')

