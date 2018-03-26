# Improve (top features) RF model based on variables analysis, class weightening
# and different sampling schemes for growing trees
# Mariana Mendoza, adapted from Ronnie Alves' codes
# February, 2012

setwd("/Users/marimendoza/Doutorado/UFRGS/Pesquisas/miRNA/RFMirTarget/")

library(gplots)
library(ellipse)
library(doParallel)

source("Scripts/functions.R")
load("Scripts/rf_model_topfeat.RData")


#-------------------------------------------------------------------------------
# SCALE DATA AND TRAIN RF MODEL

# data scaling analysis for 34 features
# boxplot -- SCALE
pdf("Plots_Feb2013/dataInspection_boxplot_NoScale.pdf")
boxplot(data.train.descr,cex.axis=0.8,cex.names=0.6,las=2)
dev.off() 

# boxplot -- NO SCALE
pdf("Plots_Feb2013/dataInspection_boxplot_Scale.pdf")
boxplot(scale(data.train.descr),cex.axis=0.8,cex.names=0.6,las=2)
dev.off()


# tune scaled RF model and perform repeated CV
data.train.descr.scaled <- scale(data.train.descr)
data.train.descr.scaled[is.nan(data.train.descr.scaled)] <- 0 

clus<- makeCluster(spec=4,type='PSOCK')
registerDoParallel(clus)
rf.scaled.grid <- expand.grid(.mtry=(1:(ncol(data.train.descr.scaled[,topFeat])))) #number of predictors (mtry) to test
rf.scaled.fit <- TuneModel(data=data.train.descr.scaled[,topFeat],labels=data.train.class,method="rf",metric="ROC",myFolds=myFolds,number=10,repeats=5,myGrid=rf.scaled.grid)
stopCluster(clus)
rm(clus)

# train model with optimized 'mtry' value (using complete data set) 
data.scaled <- data.frame(data.train.descr.scaled[,topFeat],data.train.class)
colnames(data.scaled)[dim(data.scaled)[2]] <- "Class"
rf.scaled.model <- TrainRFModel(data.scaled ,numTrees=500,numVar=rf.scaled.fit$bestTune[,'.mtry'],seed=9)

# model performance
rf.scaled.performance <- ComputeOOBStatistics(rf.scaled.model,verbose=TRUE)

# estimate model performance in terms of a confusion matrix from repeated CV
oof <- rf.scaled.fit$pred
oof <- oof[oof$.mtry==rf.scaled.fit$bestTune[,'.mtry'],]
repeats<-rf.scaled.fit$control$repeats
rf.scaled.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats)
rm(oof,repeats)
rm(data.scaled)


#---

# BUILD TOP FEATURES MODEL ELIMITATING CORRELATED VARIABLES

# analysis of correlation
pdf("Plots_Feb2013/dataInspection_correlation_NoScale.pdf")
#setEPS()
#postscript("Plots_Feb2013/dataInspection_correlation_NoScale.eps",height=9,width=9)
nzv <- nearZeroVar(data.train.descr)
data.train.descr <- data.train.descr[, -nzv]
corr.data <- cor(data.train.descr)
#ord <- order(corr.data[1,])
xc <- corr.data#[ord, ord]
colors <- c("#A50F15","#DE2D26","#FB6A4A","#FCAE91","#FEE5D9","white",
            "#EFF3FF","#BDD7E7","#6BAED6","#3182BD","#08519C")   
plotcorr(xc, col=colors[5*xc + 6],cex.lab=0.7)
dev.off()
rm(corr.data,ord,xc,colors)

# eliminate features with strong correlation (visual analysis)
# seedAU and totalGC
featuresOk <- setdiff(colnames(data.train.descr.topFeat),c("seedAU","totalGC"))

clus<- makeCluster(spec=4,type='PSOCK')
registerDoParallel(clus)
rf.nocorr.grid <- expand.grid(.mtry=(1:(ncol(data.train.descr[,featuresOk])))) #number of predictors (mtry) to test
rf.nocorr.fit <- TuneModel(data.train.descr[,featuresOk],labels=data.train.class,method="rf",metric="ROC",myFolds=myFolds,number=10,repeats=5,myGrid=rf.nocorr.grid)
stopCluster(clus)
rm(clus)


data.nocorr <- data.frame(data.train.descr[,featuresOk],data.train.class)
colnames(data.nocorr)[dim(data.nocorr)[2]] <- "Class"

# train model with optimized 'mtry' value (using complete data set) 
rf.nocorr.model <- TrainRFModel(data.nocorr,numTrees=500,numVar=rf.nocorr.fit$bestTune[,'.mtry'],seed=9)

# model performance
rf.nocorr.performance <- ComputeOOBStatistics(rf.nocorr.model,verbose=TRUE)

# estimate model performance in terms of a confusion matrix from repeated cross-validation
oof <- rf.nocorr.fit$pred
oof <- oof[oof$.mtry==rf.nocorr.fit$bestTune[,'.mtry'],]
repeats<-rf.nocorr.fit$control$repeats
rf.nocorr.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats)
rm(oof,repeats)

#---

# BUILD TOP FEATURES MODEL ELIMITATING CORRELATED VARIABLES AND APPLYING CLASS WEIGHTS

data.nocorr <- data.frame(data.train.descr[,featuresOk],data.train.class)
colnames(data.nocorr)[dim(data.nocorr)[2]] <- "Class"

clus<- makeCluster(spec=4,type='PSOCK')
registerDoParallel(clus)
rf.nocorr.wgt.grid <- expand.grid(.mtry=(1:(ncol(data.train.descr[,featuresOk])))) #number of predictors (mtry) to test
rf.nocorr.wgt.fit <- TuneModel(data.train.descr[,featuresOk],labels=data.train.class,method="rf",metric="ROC",myFolds=myFolds,number=10,repeats=5,myGrid=rf.nocorr.grid,weights=c(.70,.30))
stopCluster(clus)
rm(clus)

# train model with optimized 'mtry' value (using complete data set) 
rf.nocorr.wgt.model <- TrainRFModel(data.nocorr,numTrees=500,numVar=rf.nocorr.wgt.fit$bestTune[,'.mtry'],seed=9,weights=c(.70,.30))

# model performance
rf.nocorr.wgt.performance <- ComputeOOBStatistics(rf.nocorr.wgt.model,verbose=TRUE)

# estimate model performance in terms of a confusion matrix from repeated cross-validation
oof <- rf.nocorr.wgt.fit$pred
oof <- oof[oof$.mtry==rf.nocorr.wgt.fit$bestTune[,'.mtry'],]
repeats<-rf.nocorr.wgt.fit$control$repeats
rf.nocorr.wgt.performance.cv <- EstimatePerformanceCV(oof=oof,repeats=repeats)
rm(oof,repeats)

rm(data.nocorr)

#---

save(rf.scaled.fit,rf.scaled.grid,rf.scaled.model,rf.scaled.performance,rf.scaled.performance.cv,rf.nocorr.fit,rf.nocorr.grid,rf.nocorr.model,rf.nocorr.performance,rf.nocorr.performance.cv,rf.nocorr.wgt.fit,rf.nocorr.wgt.grid,rf.nocorr.wgt.model,rf.nocorr.wgt.performance,rf.nocorr.wgt.performance.cv,file="Scripts/rf_model_improve.RData")



