# Test trained models using independent test set
# Mariana Mendoza
# February, 2012

setwd("/Users/marimendoza/Doutorado/UFRGS/Pesquisas/miRNA/RFMirTarget/")
source("Scripts/functions.R")
source("Scripts/loadData.R")

load("Scripts/rf_model_34feat.RData")
load("Scripts/rf_model_topfeat.RData")
load("Scripts/rf_model_compare.RData")
load("Scripts/rf_model_compare_topFeat.RData")

library(caret)

# # remove variables with very small variance
# nzv <- nearZeroVar(data.train.descr)
# data.train.descr <- data.train.descr[, -nzv]
# data.test.descr <- data.test.descr[, -nzv]

# ---- LOAD TEST DATA ---- 

# # load testing data - Tarbase
data.test <- read.table("Data/testset_multimitar_rnaduplex.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE)
x<-as.character(data.test[,ncol(data.test)])
x[which(x=="Non-Target")]<-"NonTarget"
data.test[,ncol(data.test)]<-factor(x,labels=c("NonTarget","Target"))

data.test.descr <- data.test[,1:(ncol(data.test)-1)]
data.test.class <- data.test[,ncol(data.test)]

rm(x)


# ---- VARIABLES ---- 

methods.names <- c("RF","SVM","NB","KNN","GLM")#,"J48")
color <- rainbow(length(methods.names))   
models.allfeat <- list( RF=rf.fit, SVM=svm.fit, NB=nb.fit, KNN=knn.fit, GLM=glm.fit)#,J48=j48.fit)
models.top <- list( RF=rf.top.fit, SVM=svm.top.fit, NB=nb.top.fit, KNN=knn.top.fit,GLM=glm.top.fit)# ,J48=j48.top.fit)


# PREDICTION BASED ON ALL FEATURES
rf.predict.test <- predict(rf.fit,data.test.descr,type="prob")
svm.predict.test <- predict(svm.fit,data.test.descr,type="prob")
knn.predict.test <- predict(knn.fit,data.test.descr,type="prob")
#j48.predict.test <- predict(j48.fit,data.test.descr,type="prob")
nb.predict.test <- predict(nb.fit,data.test.descr,type="prob")
glm.predict.test <- predict(glm.fit,data.test.descr,type="prob")

models.allfeat.predict <- list( RF=rf.predict.test, SVM=svm.predict.test, NB=nb.predict.test, KNN=knn.predict.test,GLM=glm.predict.test)#, J48=j48.predict.test)


# PREDICTION BASED ON TOP FEATURES
rf.top.predict.test <- predict(rf.top.fit,data.test.descr[,colnames(data.train.descr.topFeat)],type="prob")
svm.top.predict.test <- predict(svm.top.fit,data.test.descr[,colnames(data.train.descr.topFeat)],type="prob")
knn.top.predict.test <- predict(knn.top.fit,data.test.descr[,colnames(data.train.descr.topFeat)],type="prob")
#j48.top.predict.test <- predict(j48.top.fit,data.test.descr[,colnames(data.train.descr.topFeat)],type="prob")
nb.top.predict.test <- predict(nb.top.fit,data.test.descr[,colnames(data.train.descr.topFeat)],type="prob")
glm.top.predict.test <- predict(glm.top.fit,data.test.descr[,colnames(data.train.descr.topFeat)],type="prob")

models.top.predict <- list( RF=rf.top.predict.test, SVM=svm.top.predict.test, NB=nb.top.predict.test, KNN=knn.top.predict.test, GLM=glm.top.predict.test)#,J48=j48.top.predict.test)


#save(models.allfeat,models.allfeat.predict,models.top,models.top.predict,methods.names,file="compare_indtestset_all.RData")


# ---- EVALUATION ---- 

# compute AUC score and plot ROC curve (if plotROC=TRUE)

auc.allFeat <- CompareModelsAUCScore(models.allfeat,data.test.descr,data.test.class,plotROC=TRUE,filename="figure_roc_testset.pdf")
auc.top <- CompareModelsAUCScore(models.top,data.test.descr[,colnames(data.train.descr.topFeat)],data.test.class,plotROC=TRUE,filename="figure_roc_testset_top.pdf")

#--- 

# compute statistics based on confusion matrix
performance.allFeat <- CompareModelsConfusionMatrix(models.allfeat,data.test.descr,data.test.class)
performance.top <- CompareModelsConfusionMatrix(models.top,data.test.descr[,colnames(data.train.descr.topFeat)],data.test.class)

#--- 

# perform permutation test and compute p-values

pvalues.allFeat <- NULL
for (jj in 1:length(models.allfeat.predict)) {
  pvalue<-PermutationTest(models.allfeat.predict[[jj]],data.test.class,auc.allFeat[jj],N=2000)
  pvalues.allFeat <- c(pvalues.allFeat,pvalue)
}
names(pvalues.allFeat) <- names(models.allfeat.predict)
rm(pvalue)

pvalues.top <- NULL
for (jj in 1:length(models.top.predict)) {
  pvalue<-PermutationTest(models.top.predict[[jj]],data.test.class,auc.top[jj],N=2000)
  pvalues.top <- c(pvalues.top,pvalue)
}
names(pvalues.top) <- names(models.top.predict)
rm(pvalue)

#--- 

# For each algorithm we observe whether the true positive predictions 
# are uniformly distributed along the entire ranked list or distributed 
# preferentially at the top of the ranked list.

o<-order(svm.top.predict.test[,"Target"],decreasing=TRUE)
data.class.ordered <- data.test.class[o]
pred<-svm.top.predict.test[o,"Target"]
pred[which(pred>0.5)] <- "Target"
pred[which(pred<=0.5)] <- "NonTarget"
pred.discrete <- rep(0,length(pred))
for (ii in 1:length(pred)) {
  if (pred[ii] == "Target" && data.class.ordered[ii] == "Target") {
    pred.discrete[ii] <- 1
  }
}
x<-(pred == data.class.ordered)
#---

We found that, 78.95% of the total true positive
predictions fall within 20th percentile of MultiMiTar ranked list;
this is very high compared to the other algorithms used in this
study. For example, only 17.65%, 30.77%, 25.93%, 25%, 36.36%
and 55% true positive predictions lie within the 20th percentile of
the ranked list of miRanda, NBmiRTar, PITA, RNAhybrid,
DIANA-microT 3.0 and TargetMiner, respectively (see Figure 4).
This clearly elucidates the fact that the ranking provided by
MultiMiTar is superior to the ranking provided by the existing
target prediction algorithms including TargetMiner



# -- 

# plot ROC curves

# add Miranda performance
points((1-0.33),0.79,pch=22) # SEN - 0.79, SPE - 0.33

# add TargetSpy-noseed-sens performance
points((1-0.54),0.45,pch=23) # SEN - 0.45, SPE - 0.54

# add TargetSpy-seed-sens performance
points((1-0.78),0.377,pch=24) # SEN - 0.377, SPE - 0.78

# add RF TOP 12 performance
points((1-0.38),0.8410,pch=25) # SEN - 0.8410, SPE - 0.38

# add SVM TOP 12 performance
points((1-0.30),0.8468,pch=21) # SEN - 0.8468, SPE - 0.30

# ----

predProbs<-extractProb(models,testX = data.test.descr,  testY = data.test.class)
predProbs.top<-extractProb(models.top, testX = data.test.descr[,colnames(data.train.descr.topFeat)],  testY = data.test.class)


#---

#save(rf.predict.test,rf.predict.test.pred,svm.predict.test,svm.predict.test.pred,knn.predict.test,knn.predict.test.pred,j48.predict.test,j48.predict.test.pred,nb.predict.test,nb.predict.test.pred,glm.predict.test,glm.predict.test.pred,models,predProbs,file="rf_compare_indeptestset.RData")

#save(rf.top.predict.test,rf.top.predict.test.pred,svm.top.predict.test,svm.top.predict.test.pred,knn.top.predict.test,knn.top.predict.test.pred,j48.top.predict.test,j48.top.predict.test.pred,nb.top.predict.test,nb.top.predict.test.pred,glm.top.predict.test,glm.top.predict.test.pred,models.top,predProbs.top,file="rf_compare_indeptestset_top.RData")
