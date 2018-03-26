# Train other classifiers with the set of 34 features
# Compare with RFMirTarget
# Mariana Mendoza
# February, 2012

setwd("/Users/marimendoza/Doutorado/UFRGS/Pesquisas/miRNA/RFMirTarget/")
source("Scripts/functions.R")
source("Scripts/loadData.R")

library(randomForest)
library(caret)
library(doParallel)


#-------------------------------------------------------------------------------
#TUNE AND TRAIN SEVERAL OTHER CLASSIFIERS AND COMPARE AGAINST RF

# configure repeated CV for train() function
myControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  repeats = 5,
  index = myFolds,
  savePredictions=TRUE,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  selectionFunction = "oneSE")


# load RF model -  repeated (5) 10 fold CV 
load("Scripts/rf_model_34feat.RData")

#-------------------------------------------------------------------------------
# 
# # remove variables with very small variance
# nzv <- nearZeroVar(data.train.descr)
# data.train.descr <- data.train.descr[, -nzv]
# data.test.descr <- data.test.descr[, -nzv]

# 34 FEATURES

# k-nearest neighbor (KNN)
clus <- makeCluster(spec=4,type='PSOCK')
registerDoParallel(clus)
knn.fit <- train(data.train.descr,data.train.class,
                 method = "knn",
                 metric="ROC",
                 trControl = myControl,
                 tuneGrid=expand.grid(expand.grid(.k=5:50)))
stopCluster(clus)
oof <- knn.fit$pred
repeats <- knn.fit$control$repeats
knn.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats)
rm(oof,repeats)


# decision trees (j48)
clus <- makeCluster(spec=4,type='PSOCK')
registerDoParallel(clus)
j48.fit <- train(data.train.descr,data.train.class,
                 method = "J48",
                 metric="ROC",
                 trControl = myControl)
stopCluster(clus)
detach("package:RWeka", unload=TRUE)
oof <- j48.fit$pred
repeats <- j48.fit$control$repeats
j48.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats)
rm(oof,repeats)


# naiveBayes
clus <- makeCluster(spec=4,type='PSOCK')
registerDoParallel(clus)
nb.fit <- train(data.train.descr,data.train.class, 
                method = "nb",
                metric="ROC",
                trControl = myControl)
stopCluster(clus)
detach("package:klaR", unload=TRUE)
oof <- nb.fit$pred
repeats <- nb.fit$control$repeats
nb.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats)
rm(oof,repeats)


# SVM
clus <- makeCluster(spec=4,type='PSOCK')
registerDoParallel(clus)
svm.fit <- train(data.train.descr,data.train.class, 
                 method = "svmRadial",
                 tuneLength = 10,
                 metric="ROC",
                 trControl = myControl, 
                 scaled = FALSE)
stopCluster(clus)
detach("package:kernlab", unload=TRUE)
oof <- svm.fit$pred
repeats <- svm.fit$control$repeats
svm.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats)
rm(oof,repeats)


# # 2 multilayer perceptron (MLP)
# clus <- makeCluster(spec=4,type='PSOCK')
# registerDoParallel(clus)
# nnet.fit <- train(data.train.descr,data.train.class,
#                   method='mlp',
#                   metric="ROC",
#                   tuneGrid=expand.grid(.size=1:10),
#                   trControl=myControl)
# stopCluster(clus)
# detach("package:RSNNS", unload=TRUE)
# oof <- nnet.fit$pred
# repeats <- nnet.fit$control$repeats
# nnet.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats)
# rm(oof,repeats)

                 
# Linear Model
clus <- makeCluster(spec=4,type='PSOCK')
registerDoParallel(clus)
glm.fit <- train(data.train.descr,data.train.class,
                  method='glm',
                  metric="ROC",
                  trControl=myControl)
stopCluster(clus)
oof <- glm.fit$pred
repeats <- glm.fit$control$repeats
glm.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats)
rm(oof,repeats)

#---

#save(knn.fit,knn.grid,knn.performance.cv,j48.fit,j48.performance.cv,nb.fit,nb.performance.cv,svm.fit,svm.performance.cv, nnet.fit,nnet.performance.cv,glm.fit,glm.performance.cv,myControl,file="Scripts/rf_model_compare.RData")
                 
#rm(knn.fit,knn.grid,knn.performance.cv,j48.fit,j48.performance.cv,nb.fit,nb.performance.cv,svm.fit,svm.performance.cv, nnet.fit,nnet.performance.cv,glm.fit,glm.performance.cv,myControl)

#-------------------------------------------------------------------------------
                 
# TOP FEATURES

# k-nearest neighbor (KNN)
clus <- makeCluster(spec=4,type='PSOCK')
registerDoParallel(clus)
knn.top.fit <- train(data.train.descr.topFeat,data.train.class,
                 method = "knn",
                 metric="ROC",
                 trControl = myControl,
                 tuneGrid=expand.grid(expand.grid(.k=5:50)))
stopCluster(clus)
oof <- knn.top.fit$pred
repeats <- knn.top.fit$control$repeats
knn.top.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats)
rm(oof,repeats)


# decision trees (j48)
clus <- makeCluster(spec=4,type='PSOCK')
registerDoParallel(clus)
j48.top.fit <- train(data.train.descr.topFeat,data.train.class,
                 method = "J48",
                 metric="ROC",
                 trControl = myControl)
stopCluster(clus)
detach("package:RWeka", unload=TRUE)
oof <- j48.top.fit$pred
repeats <- j48.top.fit$control$repeats
j48.top.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats)
rm(oof,repeats)



# naiveBayes
clus <- makeCluster(spec=4,type='PSOCK')
registerDoParallel(clus)
nb.top.fit <- train(data.train.descr.topFeat,data.train.class, 
                method = "nb",
                metric="ROC",
                trControl = myControl)
stopCluster(clus)
detach("package:klaR", unload=TRUE)
oof <- nb.top.fit$pred
repeats <- nb.top.fit$control$repeats
nb.top.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats)
rm(oof,repeats)




# SVM
clus <- makeCluster(spec=4,type='PSOCK')
registerDoParallel(clus)
svm.top.fit <- train(data.train.descr.topFeat,data.train.class, 
                 method = "svmRadial",
                 tuneLength = 10,
                 metric="ROC",
                 trControl = myControl, 
                 scaled = FALSE)
stopCluster(clus)
detach("package:kernlab", unload=TRUE)
oof <- svm.top.fit$pred
repeats <- svm.top.fit$control$repeats
svm.top.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats)
rm(oof,repeats)



# # 2 multilayer perceptron (MLP)
# clus <- makeCluster(spec=4,type='PSOCK')
# registerDoParallel(clus)
# nnet.top.fit <- train(data.train.descr.topFeat,data.train.class,
#                   method='mlp',
#                   metric="ROC",
#                   tuneGrid=expand.grid(.size=1:10),
#                   trControl=myControl)
# stopCluster(clus)
# detach("package:RSNNS", unload=TRUE)
# oof <- nnet.top.fit$pred
# repeats <- nnet.top.fit$control$repeats
# nnet.top.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats)
# rm(oof,repeats)

                 
# Linear Model
clus <- makeCluster(spec=4,type='PSOCK')
registerDoParallel(clus)
glm.top.fit <- train(data.train.descr.topFeat,data.train.class,
                 method='glm',
                 metric="ROC",
                 trControl=myControl)
stopCluster(clus)
oof <- glm.top.fit$pred
repeats <- glm.top.fit$control$repeats
glm.top.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats)
rm(oof,repeats)
                 
#---
                
#save(knn.top.fit,knn.top.performance.cv,j48.top.fit,j48.top.performance.cv,nb.top.fit,nb.top.performance.cv,svm.top.fit,svm.top.performance.cv,glm.top.fit,glm.top.performance.cv,myControl,file="Scripts/rf_model_compare_topFeat.RData")
                 
#rm(knn.top.fit,knn.top.performance.cv,j48.top.fit,j48.top.performance.cv,nb.top.fit,nb.top.performance.cv,svm.top.fit,svm.top.performance.cv, nnet.top.fit,nnet.top.performance.cv,glm.top.fit,glm.top.performance.cv,myControl)