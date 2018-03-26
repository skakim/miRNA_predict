# Compare classifiers in terms of confidence intervals (95%) over AUC scores
# Mariana Mendoza
# February, 2012

#setwd("/Users/marimendoza/Doutorado/UFRGS/Pesquisas/miRNA/RFMirTarget/")
source("Scripts/functions.R")
load("Scripts/rf_model_34feat.RData")
load("Scripts/rf_model_topfeat.RData")
load("Scripts/rf_model_compare.RData")

library(caret)

# compare models all features
resamps <- resamples(list(rf=rf.fit, knn=knn.fit, svm=svm.fit,nb=nb.fit,j48=j48.fit, glm=glm.fit))
methods <- resamps$models
metric <- "ROC"

resamps.val <- resamps$values

numMethods <- length(methods)
numSamples <- dim(resamps.val)[1]
index <- paste(methods,"~",metric,sep="")

all.means <- apply(resamps.val[index], 2, mean)
all.sd <- apply(resamps.val[index], 2, sd)

# order results to plot in decreasing format
o<-order(all.means,decreasing=TRUE)

pdf("Plots_Feb2013/Figure_auc_barplot.pdf")
barx <- barplot2(all.means[o], beside=TRUE,names.arg=methods[o],ylim=c(0,1.0), col="blue", axis.lty=1, xlab="Method", ylab="Average AUC score")
error.bar(barx,all.means[o], all.sd[o])
box()
# text(c, all.means[o]+all.sd[o]+0.0001, labels = format(all.means[o], digits=3),pos = 3, cex = .85)
dev.off()




# compare models top features
resamps <- resamples(list(rf=rf.top.fit, knn=knn.top.fit, svm=svm.top.fit,nb=nb.top.fit,j48=j48.top.fit, glm=glm.top.fit))
methods <- resamps$models
metric <- "ROC"

resamps.val <- resamps$values

numMethods <- length(methods)
numSamples <- dim(resamps.val)[1]
index <- paste(methods,"~",metric,sep="")

all.means <- apply(resamps.val[index], 2, mean)
all.sd <- apply(resamps.val[index], 2, sd)

# order results to show in decreasing format
o<-order(all.means,decreasing=TRUE)


pdf("Plots_Feb2013/Figure_auc_barplot_top.pdf")
barx <- barplot2(all.means[o], beside=TRUE,names.arg=methods[o],ylim=c(0,1.0), col="blue", axis.lty=1, xlab="Method", ylab="Average AUC score")

error.bar(barx,all.means[o], all.sd[o])
box()
dev.off()


# ---- temporary code
# resamps <- resamples(list(rf=rf.fit, knn=knn.fit,glm=glm.fit,svm=svm.fit,nb=nb.fit,j48=j48.fit))#, nnet=fit.nnet))
# print(summary(resamps))
# 
# # compute estimates of difference and p-values for pairwise comparison among classifiers
# diffs <- diff(resamps,metric="ROC")
# print(summary(diffs))
# 
# 
# bwplot(resamps, metric="ROC")
# dotplot(diffs)
# dotplot(resamps,metric="ROC")
# pdf("Plots_Feb2013/Figure_auc_diffs_top.pdf")
# diffs <- diff(resamps,metric="ROC")
# levelplot(diffs, what = "differences",region=TRUE)
# dev.off()

# ---- temporary code

