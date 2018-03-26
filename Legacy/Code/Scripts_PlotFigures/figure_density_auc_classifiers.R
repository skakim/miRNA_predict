# Compare classifiers in terms of density curves
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
pdf("Plots_Feb2013/Figure_auc_density.pdf")
PlotDensityCurves(methods,"ROC",resamps$values,rug=TRUE)
dev.off()


# compare models top features
resamps.top <- resamples(list(rf=rf.top.fit, knn=knn.top.fit, svm=svm.top.fit,nb=nb.top.fit,j48=j48.top.fit, glm=glm.top.fit))
methods <- resamps.top$models
pdf("Plots_Feb2013/Figure_auc_density_top.pdf")
PlotDensityCurves(methods,"ROC",resamps.top$values,rug=TRUE)
dev.off()
