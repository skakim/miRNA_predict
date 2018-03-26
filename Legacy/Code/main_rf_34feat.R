# Train RF model with the set of 34 features.
# Test using repeated 10-fold CV
# Mariana Mendoza
# December, 2012

#setwd("/Users/marimendoza/Doutorado/UFRGS/Pesquisas/miRNA/RFMirTarget/")
#source("Scripts/functions.R")
#source("Scripts/loadData.R")

library(randomForest)
library(caret)
library(doParallel)


#-------------------------------------------------------------------------------
#TUNE AND TRAIN RF MODEL WITH 34 FEATURES

# remove variables with very small variance
nzv <- nearZeroVar(data.train.descr)
data.train.descr <- data.train.descr[, -nzv]
data.test.descr <- data.test.descr[, -nzv]

# tune mtry (number of predictors) and perform cross-validation on data
# use doParallel package to run in parallel 
clus<- makeCluster(spec=4,type='PSOCK')
registerDoParallel(clus)
rf.grid <- expand.grid(.mtry=(5:(ncol(data.train.descr)))) #number of predictors (mtry) to test
rf.fit <- TuneModel(data=data.train.descr,labels=data.train.class,method="rf",metric="ROC",myFolds=myFolds,number=10,repeats=5,myGrid=rf.grid)
optimal.mtry <- rf.fit$bestTune[,'.mtry']
stopCluster(clus)
rm(clus)

# train model with optimized 'mtry' value (using complete data set) 
rf.model <- TrainRFModel(data.train,numTrees=500,numVar=optimal.mtry,seed=9)

# model performance
rf.performance <- ComputeOOBStatistics(rf.model,verbose=TRUE)

# features importance
rf.featImp <- ComputeFeaturesImportance(rf.model, criteria="gini")
rf.performance.rfs <- RestrictedForwardSelection(data.train,rf.model,rf.featImp,verbose=TRUE)

# estimate model performance in terms of a confusion matrix from repeated cross-validation
oof <- rf.fit$pred
oof <- oof[oof$.mtry==rf.fit$bestTune[,'.mtry'],]
repeats<-rf.fit$control$repeats
rf.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats)

# estimate error rates from repeated cross-validation
oof <- rf.fit$pred
oof <- oof[oof$.mtry==rf.fit$bestTune[,'.mtry'],]
repeats<-rf.fit$control$repeats
rf.error.cv <- EstimateErrorRateCV(oof=oof,repeats=repeats)
rm(oof,repeats)

save(rf.featImp,rf.performance,rf.grid,rf.model,rf.performance.cv,rf.fit,rf.performance.rfs,rf.predict.test,optimal.mtry,nzv,file="Scripts/rf_model_34feat.RData")