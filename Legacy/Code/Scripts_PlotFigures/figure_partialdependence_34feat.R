# Compare classifiers in terms of confidence intervals (95%) over AUC scores
# Mariana Mendoza
# February, 2012

#setwd("/Users/marimendoza/Doutorado/UFRGS/Pesquisas/miRNA/RFMirTarget/")
#source("Scripts/functions.R")
load("Scripts/rf_model_34feat.RData")
load("Scripts/rf_model_topfeat.RData")

library(caret)
library(randomForest)

pdf("partialDependence_34featmodel.pdf")
op <- par(mfrow=c(4, 3))
# plot the partial dependence: marginal effect of a variable on the class probability (classification) 
for (ii in 1:length(topFeat)) {
  partialPlot(rf.model,data.train,colnames(data.train)[topFeat[ii]],"Target",xlab=colnames(data.train)[topFeat[ii]],ylab = "Partial Dependence",main="")
}
par(op)
dev.off()

pdf("partialDependence_12featmodel.pdf")
op <- par(mfrow=c(4, 3))
# plot the partial dependence: marginal effect of a variable on the class probability (classification) 
for (ii in 1:length(topFeat)) {
  partialPlot(rf.top.model,data.train,colnames(data.train)[topFeat[ii]],"Target",xlab=colnames(data.train)[topFeat[ii]],ylab = "Partial Dependence",main="")
}
par(op)
dev.off()



