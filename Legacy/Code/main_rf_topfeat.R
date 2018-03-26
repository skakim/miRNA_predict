# Train a RF classifier with the top features based on features importance ranking
# Mariana Mendoza
# February, 2012

#setwd("/Users/marimendoza/Doutorado/UFRGS/Pesquisas/miRNA/RFMirTarget/")
#source("Scripts/functions.R")

library(randomForest)
library(caret)
library(doParallel)


#-------------------------------------------------------------------------------
#BUILD MODEL WITH TOP RELEVANT FEATURES

# load RF model -  repeated (5) 10 fold CV 
load("Scripts/rf_model_34feat.RData")

# find maximum MCC
idx <- which(rf.performance.rfs$values[,4] == max(rf.performance.rfs$values[,4]))
topFeat <- rf.performance.rfs$features[1:idx]
data.train.descr.topFeat <- data.train.descr[,topFeat]

clus<- makeCluster(spec=4,type='PSOCK')
registerDoParallel(clus)
rf.top.grid <- expand.grid(.mtry=(1:(ncol(data.train.descr.topFeat)))) #number of predictors (mtry) to test
rf.top.fit <- TuneModel(data=data.train.descr.topFeat,labels=data.train.class,method="rf",metric="ROC",myFolds=myFolds,number=10,repeats=5,myGrid=rf.top.grid)
stopCluster(clus)
rm(clus)

# train model with optimized 'mtry' value (using complete data set) 
rf.top.model <- TrainRFModel(data.train[,c(topFeat,35)],numTrees=500,numVar=rf.top.fit$bestTune[,'.mtry'],seed=9)

# model performance
rf.top.performance <- ComputeOOBStatistics(rf.model,verbose=TRUE)

# estimate model performance (top 12 features) in terms of a confusion matrix from repeated cross-validation
oof <- rf.top.fit$pred
oof <- oof[oof$.mtry==rf.top.fit$bestTune[,'.mtry'],]
repeats<-rf.top.fit$control$repeats
rf.top.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats)
rm(oof,repeats)

# estimate error rates from repeated cross-validation
oof <- rf.top.fit$pred
oof <- oof[oof$.mtry==rf.top.fit$bestTune[,'.mtry'],]
repeats <- rf.top.fit$control$repeats
rf.top.error.cv <- EstimateErrorRateCV(oof=oof,repeats=repeats)
rm(oof,repeats)

#---

save(topFeat,rf.top.fit,rf.top.grid,rf.top.model,rf.top.performance,rf.top.performance.cv,data.train.descr.topFeat,file="Scripts/rf_model_topfeat.RData")

