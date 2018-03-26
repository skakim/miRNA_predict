# Compare classifiers in terms of confidence intervals (95%) over AUC scores
# Mariana Mendoza
# February, 2012

#setwd("/Users/marimendoza/Doutorado/UFRGS/Pesquisas/miRNA/RFMirTarget/")
#source("Scripts/functions.R")
load("Scripts/main_rf_34feat.R")
load("Scripts/main_rf_topfeat.R")
load("Scripts/main_compare_classifiers.R")

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

o<-order(all.means,decreasing=TRUE)

ci.up <- all.means[o] + (1.96*(all.sd[o]/sqrt(numSamples)))
ci.low <- all.means[o] - (1.96*(all.sd[o]/sqrt(numSamples)))

pdf("Plots_Feb2013/Figure_auc_ci.pdf")
plotCI(all.means[o],ui=ci.up,li=ci.low,gap=0.3,lwd=2,sfrac=0.005,xlab="Methods",ylab="AUC (95% confidence intervals)",horizontal=TRUE)
dev.off()

pdf("Plots_Feb2013/Figure_auc_ci2.pdf")
dotplot(resamps,metric="ROC")
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

o<-order(all.means,decreasing=TRUE)

ci.up <- all.means[o] + (1.96*(all.sd[o]/sqrt(numSamples)))
ci.low <- all.means[o] - (1.96*(all.sd[o]/sqrt(numSamples)))

pdf("Plots_Feb2013/Figure_auc_ci_top.pdf")
plotCI(all.means[o],ui=ci.up,li=ci.low,gap=0.3,lwd=2,sfrac=0.005,xlab="Methods",ylab="AUC (95% confidence intervals)")
dev.off()

pdf("Plots_Feb2013/Figure_auc_ci_top2.pdf")
dotplot(resamps,metric="ROC")
dev.off()


