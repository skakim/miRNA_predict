
library("plyr")

setwd("/home/I864741/TCC/Analysis/")

##########################################################
##    create a mapping function to standardize genes    ##
##                names using HGCN data                 ##
##########################################################

####  Create mapping structure  to standardize genes names ####
hgcn.nomenclature <-read.delim("/home/I864741/TCC/Analysis/IDs_and_interaction_data/HGCN_completeSet.txt",stringsAsFactors = FALSE)

genename.map <- NULL
names <- NULL
for (ii in 1:dim(hgcn.nomenclature)[1]) {
    ## treat alias_symbol field
    temp <- hgcn.nomenclature[ii,"alias_symbol"]
    if (temp != "") {
        x <- unlist(strsplit(temp,"\\|"))
        names <- c(names,x)
        genename.map <- c(genename.map ,rep(hgcn.nomenclature[ii,"symbol"],length(x)))
    }
    names <- c(names,hgcn.nomenclature[ii,"symbol"])
    genename.map <- c(genename.map ,hgcn.nomenclature[ii,"symbol"])
    
    ## treat prev_symbol field
    temp <- hgcn.nomenclature[ii,"prev_symbol"]
    if (temp != "") {
        x <- unlist(strsplit(temp,"\\|"))
        names <- c(names,x)
        genename.map <- c(genename.map ,rep(hgcn.nomenclature[ii,"symbol"],length(x)))
    }
}
names(genename.map) <- names



## include symbols manually fixed
manualfix <- read.csv("/home/I864741/TCC/Analysis/IDs_and_interaction_data/HGCN_symbols_included.csv",stringsAsFactors = FALSE,header=FALSE)

names <- manualfix[,1]
temp <- manualfix[,2]
names(temp) <- names
genename.map <- c(genename.map,temp)


## Mapping function
standardizeGeneSymbols <- function(geneList,mapFunction,uppercase=FALSE){
    if(uppercase) {
        geneList <- toupper(geneList)
    }
    temp<-mapFunction[as.character(geneList)]
    temp <- as.character(temp)
    return(temp)
}
rm(names,temp,x,hgcn.nomenclature,ii,manualfix)

#save(genename.map,file="~/Dropbox/UFRGS/Orientandos/TCC/2018_GabrielMoita_ensembleMirnas/HGCN_MapFunction.RData")


##########################################################
##    create a mapping function to standardize miRNA    ##
##               names using miRBase v21                ##
##########################################################

##### miRNA mapping structure ######

## load data
mirna.mature.v21 <- read.delim("/home/I864741/TCC/Analysis/IDs_and_interaction_data/mirna_mature_v21.txt", header=FALSE)
mirna.mature.v21 <- mirna.mature.v21[grep("hsa",mirna.mature.v21[,2]),]
mirna.mature.v21[] <- apply(mirna.mature.v21,2,as.character)

## create mapping function
map.mirnas.v21 <-NULL
names.prev <- NULL
for(ii in 1:nrow(mirna.mature.v21)) {
    temp <- unlist(strsplit(mirna.mature.v21[ii,3],";"))
    names.prev <- c(names.prev,temp)
    names.prev <- c(names.prev,mirna.mature.v21[ii,2])
    map.mirnas.v21 <- c(map.mirnas.v21,rep(mirna.mature.v21[ii,2],(length(temp)+1)))
    rm(temp)
}
names(map.mirnas.v21) <- names.prev
rm(names.prev)


#save(map.mirnas.v21,file="~/Dropbox/UFRGS/Orientandos/TCC/2018_GabrielMoita_ensembleMirnas/miRBasev21_MapFunction.RData")


load("/home/I864741/TCC/Analysis/processInteractionData_output/miRBasev21_MapFunction.RData")
load("/home/I864741/TCC/Analysis/processInteractionData_output/HGCN_MapFunction.RData")

## Mapping function
standardizeGeneSymbols <- function(geneList,mapFunction,uppercase=FALSE){
    if(uppercase) {
        geneList <- toupper(geneList)
    }
    temp<-mapFunction[as.character(geneList)]
    temp <- as.character(temp)
    return(temp)
}


##### TarBase v7.0 ##### 

## load and pre-process data
tarbase <- read.csv("/home/I864741/TCC/Analysis/IDs_and_interaction_data/TarBase_v7.0.csv",header=TRUE,stringsAsFactors = FALSE,sep="\t")
tarbase <- tarbase[which(tarbase$species=="Homo sapiens"),]

## map miRNAs name according to miRBase v21 
tarbase$mirna <- as.character(map.mirnas.v21[as.character(tarbase$mirna)])

## map gene symbols to HGCN standard
tarbase$geneName <- standardizeGeneSymbols(tarbase$geneName,genename.map,uppercase = FALSE)

## remove rows with NAs (could not be correctly mapped)
tarbase <- tarbase[-union(which(is.na(tarbase$geneName)),which(is.na(tarbase$mirna))),]


## separate positive examples (POSITIVE, DIRECT, DOWN)
tarbase.positive <- tarbase[which(tarbase$positive_negative=="POSITIVE"),]
tarbase.positive <- tarbase.positive[which(tarbase.positive$direct_indirect=="DIRECT"),]
tarbase.positive <- tarbase.positive[which(tarbase.positive$up_down=="DOWN"),]

## remove HITS-CLIP and PAR-CLIP: are considered "limited evidence"
tarbase.positive<- tarbase.positive[-which(tarbase.positive$method %in% c("HITS-CLIP","PAR-CLIP")),]

#remove duplicates of mirna-gene pairs, and filter columns
tarbase.positive <- unique(tarbase.positive[,c("mirna","geneName","geneId")])


## separate negative examples 
tarbase.negative <- tarbase[which(tarbase$positive_negative=="NEGATIVE"),]

#remove duplicates of mirna-gene pairs, and filter columns
tarbase.negative<- unique(tarbase.negative[,c("mirna","geneName","geneId")])

# tarbase.negative <- tarbase.negative[-which(tarbase.negative$method %in% c("HITS-CLIP","PAR-CLIP")),] #zero cases


# any common elements between positive and negative sets? Yes!
intersect(paste0(tarbase.positive$geneId,tarbase.positive$mirna),paste0(tarbase.negative$geneId,tarbase.negative$mirna))

# pairs of miRNA-targets that are both positive and negative examples are removed
remove.from.negatives <- which(is.element(paste0(tarbase.negative$geneId,tarbase.negative$mirna),paste0(tarbase.positive$geneId,tarbase.positive$mirna)))
remove.from.positives <- which(is.element(paste0(tarbase.positive$geneId,tarbase.positive$mirna),paste0(tarbase.negative$geneId,tarbase.negative$mirna)))
tarbase.positive <- tarbase.positive[-remove.from.positives,]
tarbase.negative <- tarbase.negative[-remove.from.negatives,]


rm(remove.from.negatives,remove.from.positives)





## Data from Release 6.1
mirtarbase <- read.delim("/home/I864741/TCC/Analysis/IDs_and_interaction_data/mirtarbase_hsa_MTI_release6.1.txt",sep="\t", header=TRUE)
mirtarbase <- mirtarbase[which(mirtarbase$Species..miRNA.=="Homo sapiens"),]


## map miRNAs name according to miRBase v21 
mirtarbase$miRNA <- as.character(map.mirnas.v21[as.character(mirtarbase$miRNA)])

## map gene symbols to HGCN standard
mirtarbase$Target.Gene <- standardizeGeneSymbols(mirtarbase$Target.Gene,genename.map,uppercase = FALSE)

## remove rows with NAs (could not be correctly mapped)
mirtarbase <- mirtarbase[-union(which(is.na(mirtarbase$miRNA)),which(is.na(mirtarbase$Target.Gene))),]

## separate positive and negative examples  (functional/non-functional)
mirtarbase.positive <- mirtarbase[which(mirtarbase$Support.Type=="Functional MTI"),]
mirtarbase.negative <- mirtarbase[which(mirtarbase$Support.Type=="Non-Functional MTI"),]

# any common elements between positive and negative sets?
intersect(paste0(mirtarbase.positive$miRNA,mirtarbase.positive$Target.Gene),paste0(mirtarbase.negative$miRNA,mirtarbase.negative$Target.Gene))


# pairs of miRNA-targets that are both positive and negative examples are removed
remove.from.negatives <- which(is.element(paste0(mirtarbase.negative$miRNA,mirtarbase.negative$Target.Gene),paste0(mirtarbase.positive$miRNA,mirtarbase.positive$Target.Gene)))
remove.from.positives <- which(is.element(paste0(mirtarbase.positive$miRNA,mirtarbase.positive$Target.Gene),paste0(mirtarbase.negative$miRNA,mirtarbase.negative$Target.Gene)))
mirtarbase.positive <- mirtarbase.positive[-remove.from.positives,]
mirtarbase.negative <- mirtarbase.negative[-remove.from.negatives,]

#remove duplicates of mirna-gene pairs, and filter columns
mirtarbase.negative <- unique(mirtarbase.negative[,c("miRNA","Target.Gene","Target.Gene..Entrez.Gene.ID.")])
mirtarbase.positive <- unique(mirtarbase.positive[,c("miRNA","Target.Gene","Target.Gene..Entrez.Gene.ID.")])


## delete original data to save space
rm(tarbase,mirtarbase)


## save workspace
save.image("preprocessed_data_mirTarBasev6.1_TarBasev7.RData")
