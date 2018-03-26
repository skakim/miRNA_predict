

load("Scripts/source_data.RData")
# 
# # load training data
# data.train <- read.table("data/trainingset_noduplicates_rnaduplex.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill=TRUE)
# x<-as.character(data.train[,ncol(data.train)])
# x[which(x=="Non-Target")]<-"NonTarget"
# data.train[,ncol(data.train)]<-factor(x,labels=c("NonTarget","Target"))
# 
# data.train.descr <- data.train[,1:(ncol(data.train)-1)]
# data.train.class <- data.train[,ncol(data.train)]
# # x<-as.character(data.train.class)
# # x[which(x=="Non-Target")]<-"NonTarget"
# # data.train.class<-factor(x,labels=c("NonTarget","Target"))
# 
# 
# # load testing data - Tarbase
# data.test <- read.table("Data/testset_tarbase_rnaduplex.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE)
# x<-as.character(data.test[,ncol(data.test)])
# x[which(x=="Non-Target")]<-"NonTarget"
# data.test[,ncol(data.test)]<-factor(x,labels=c("NonTarget","Target"))
# 
# data.test.descr <- data.test[,1:(ncol(data.test)-1)]
# data.test.class <- data.test[,ncol(data.test)]
# # x<-as.character(data.test.class)
# # x[which(x=="Non-Target")]<-"NonTarget"
# # data.test.class<-factor(x,labels=c("NonTarget","Target"))
# 
# rm(x)
# 
# # create folds for CV
# set.seed(2)
# myFolds <- createMultiFolds(y=data.train.class,k=10,times=5) #use the same fold for all classifiers
# 
# save(data.train,data.train.class,data.train.descr,data.test,data.test.class,data.test.descr,myFolds,file="Scripts/source_data.RData")