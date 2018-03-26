# Train a RF model for each of the 5 features categories
# Mariana Mendoza
# December, 2012

#setwd("/Users/marimendoza/Doutorado/UFRGS/Pesquisas/miRNA/RFMirTarget/")
#source("Scripts/functions.R")
#source("Scripts/loadData.R")

library(randomForest)
library(caret)


#-------------------------------------------------------------------------------
#TUNE AND TRAIN RF MODEL FOR EACH OF THE FEATURES GROUPS

# separe features in semantic groups
g1 <- c(1:2)      #alignment
g2 <- c(3)        #thermodynamics
g3 <- c(5:9)      #structural
g4 <- c(4,10:14)  #seed
g5 <- c(15:34)    #position-based

featNum <- length(colnames(data.train))
featCat <- list(g1,g2,g3,g4,g5)
names(featCat) <- c("alignment","thermodynamics","structural","seed","position")

#-------------------------------------------------------------------------------

# alignment features
cat <- 1
cat.size <- length(featCat[[cat]])
rf.g1.grid <- expand.grid(.mtry=(1:cat.size)) #number of predictors (mtry) to test
rf.g1.fit <- TuneModel(data=data.train.descr[,featCat[[cat]]],labels=data.train.class,method="rf",metric="ROC",myFolds=myFolds,number=10,repeats=5,myGrid=rf.g1.grid)
optimal.mtry <- rf.g1.fit$bestTune[,'.mtry']

# train model with optimized 'mtry' value (using complete data set)
rf.g1.model <- TrainRFModel(data.train[,c(featCat[[cat]],featNum)],numTrees=500,numVar=optimal.mtry,seed=9)

# model performance
rf.g1.performance <- ComputeOOBStatistics(rf.g1.model,verbose=TRUE)

oof <- rf.g1.fit$pred
oof <- oof[oof$.mtry==rf.g1.fit$bestTune[,'.mtry'],]
repeats <- rf.g1.fit$control$repeats
rf.g1.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats)
rm(oof,repeats)

#---

# thermodynamics features
cat <- 2
cat.size <- length(featCat[[cat]])
rf.g2.grid <- expand.grid(.mtry=(1:cat.size)) #number of predictors (mtry) to test
rf.g2.fit <- TuneModel(data=data.train.descr[,featCat[[cat]]],labels=data.train.class,method="rf",metric="ROC",myFolds=myFolds,number=10,repeats=5,myGrid=rf.g2.grid)
optimal.mtry <- rf.g2.fit$bestTune[,'.mtry']

# train model with optimized 'mtry' value (using complete data set)
rf.g2.model <- TrainRFModel(data.train[,c(featCat[[cat]],featNum)],numTrees=500,numVar=optimal.mtry,seed=9)

# model performance
rf.g2.performance <- ComputeOOBStatistics(rf.g2.model,verbose=TRUE)

oof <- rf.g2.fit$pred
oof <- oof[oof$.mtry==rf.g2.fit$bestTune[,'.mtry'],]
repeats <- rf.g2.fit$control$repeats
rf.g2.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats)
rm(oof,repeats)

#---

# structural features
cat <- 3
cat.size <- length(featCat[[cat]])
rf.g3.grid <- expand.grid(.mtry=(1:cat.size)) #number of predictors (mtry) to test
rf.g3.fit <- TuneModel(data=data.train.descr[,featCat[[cat]]],labels=data.train.class,method="rf",metric="ROC",myFolds=myFolds,number=10,repeats=5,myGrid=rf.g3.grid)
optimal.mtry <- rf.g3.fit$bestTune[,'.mtry']

# train model with optimized 'mtry' value (using complete data set)
rf.g3.model <- TrainRFModel(data.train[,c(featCat[[cat]],featNum)],numTrees=500,numVar=optimal.mtry,seed=9)

# model performance
rf.g3.performance <- ComputeOOBStatistics(rf.g3.model,verbose=TRUE)

oof <- rf.g3.fit$pred
oof <- oof[oof$.mtry==rf.g3.fit$bestTune[,'.mtry'],]
repeats <- rf.g3.fit$control$repeats
rf.g3.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats)
rm(oof,repeats)

#---

# seed features
cat <- 4
cat.size <- length(featCat[[cat]])
rf.g4.grid <- expand.grid(.mtry=(1:cat.size)) #number of predictors (mtry) to test
rf.g4.fit <- TuneModel(data=data.train.descr[,featCat[[cat]]],labels=data.train.class,method="rf",metric="ROC",myFolds=myFolds,number=10,repeats=5,myGrid=rf.g4.grid)
optimal.mtry <- rf.g4.fit$bestTune[,'.mtry']

# train model with optimized 'mtry' value (using complete data set)
rf.g4.model <- TrainRFModel(data.train[,c(featCat[[cat]],featNum)],numTrees=500,numVar=optimal.mtry,seed=9)

# model performance
rf.g4.performance <- ComputeOOBStatistics(rf.g4.model,verbose=TRUE)

oof <- rf.g4.fit$pred
oof <- oof[oof$.mtry==rf.g4.fit$bestTune[,'.mtry'],]
repeats <- rf.g4.fit$control$repeats
rf.g4.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats)
rm(oof,repeats)

#---

# position-based features
cat <- 5
cat.size <- length(featCat[[cat]])
rf.g5.grid <- expand.grid(.mtry=(1:cat.size)) #number of predictors (mtry) to test
rf.g5.fit <- TuneModel(data=data.train.descr[,featCat[[cat]]],labels=data.train.class,method="rf",metric="ROC",myFolds=myFolds,number=10,repeats=5,myGrid=rf.g5.grid)
optimal.mtry <- rf.g5.fit$bestTune[,'.mtry']

# train model with optimized 'mtry' value (using complete data set)
rf.g5.model <- TrainRFModel(data.train[,c(featCat[[cat]],featNum)],numTrees=500,numVar=optimal.mtry,seed=9)

# model performance
rf.g5.performance <- ComputeOOBStatistics(rf.g5.model,verbose=TRUE)

oof <- rf.g5.fit$pred
oof <- oof[oof$.mtry==rf.g5.fit$bestTune[,'.mtry'],]
repeats <- rf.g5.fit$control$repeats
rf.g5.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats)
rm(oof,repeats)

rm(cat,cat.size,optimal.mtry)

#---

save(rf.g1.fit,rf.g1.grid,rf.g1.model,rf.g1.performance,rf.g1.performance.cv,rf.g2.fit,rf.g2.grid,rf.g2.model,rf.g2.performance,rf.g2.performance.cv,rf.g3.fit,rf.g3.fit,rf.g3.grid,rf.g3.model,rf.g3.performance,rf.g3.performance.cv,rf.g4.fit,rf.g4.grid,rf.g4.model,rf.g4.performance,rf.g4.performance.cv,rf.g5.fit,rf.g5.grid,rf.g5.model,rf.g5.performance,rf.g5.performance.cv,featCat,featNum,g1,g2,g3,g4,g5,file="Scripts/rf_model_featcat.RData")

rm(rf.g1.fit,rf.g1.grid,rf.g1.model,rf.g1.performance,rf.g1.performance.cv,rf.g2.fit,rf.g2.grid,rf.g2.model,rf.g2.performance,rf.g2.performance.cv,rf.g3.fit,rf.g3.grid,rf.g3.model,rf.g3.performance,rf.g3.performance.cv,rf.g4.fit,rf.g4.grid,rf.g4.model,rf.g4.performance,rf.g4.performance.cv,rf.g5.fit,rf.g5.grid,rf.g5.model,rf.g5.performance,rf.g5.performance.cv,featCat,featNum,g1,g2,g3,g4,g5)

#---

# comparing performance of different feature categories

resamps <- resamples(list(g1=rf.g1.fit,g2=rf.g2.fit,g3=rf.g3.fit,g4=rf.g4.fit,g5=rf.g5.fit))#, nnet=fit.nnet))
print(summary(resamps))

# compute estimates of difference and p-values for pairwise comparison among classifiers
diffs <- diff(resamps,metric="ROC")
print(summary(diffs))

# use package Cairo to generate file more compatible with Inkscape
CairoFonts(
  regular="Arial:style=Regular",
  bold="Arial:style=Bold",
  italic="Arial:style=Italic",
  bolditalic="Arial:style=Bold Italic,BoldItalic",
  symbol="Symbol"
)

CairoPDF("Plots_Feb2013/Figure_scatterplot_auc_featcats.pdf",bg="white",fonts=regular)
splom(resamps, metric="ROC",col="red")


dev.off()

# parallelplot(resamps, metric="ROC")
# densityplot(resamps, metric="ROC")
# xyplot(resamps, metric="ROC")

rm(diffs,resamps)
