# Compare classifiers in terms of confidence intervals (95%) over AUC scores
# Mariana Mendoza
# February, 2012

#setwd("/Users/marimendoza/Doutorado/UFRGS/Pesquisas/miRNA/RFMirTarget/")
source("Scripts/functions.R")
load("Scripts/rf_model_34feat.RData")
load("Scripts/rf_model_topfeat.RData")
load("Scripts/rf_model_compare.RData")

library(caret)
library(gplots)

# compare models all features
resamps <- resamples(list(rf=rf.fit, knn=knn.fit, svm=svm.fit,nb=nb.fit,j48=j48.fit, glm=glm.fit))
methods <- resamps$models
metric <- "ROC"

print(summary(resamps))

# compute estimates of difference and p-values for pairwise comparison among classifiers
diffs <- diff(resamps,metric=metric)
diffs.table <- summary(diffs)$table[[1]]


pdf("Plots_Feb2013/Figure_auc_diffs.pdf")
levelplot(diffs, what = "differences",col.regions = rev(heat.colors(20)))
dev.off()

pdf("Plots_Feb2013/Figure_auc_diffs_table.pdf")
textplot(diffs.table)
dev.off()



# compare models top features
resamps <- resamples(list(rf=rf.top.fit, knn=knn.top.fit, svm=svm.top.fit,nb=nb.top.fit,j48=j48.top.fit, glm=glm.top.fit))
methods <- resamps$models
metric <- "ROC"print(summary(resamps))

# compute estimates of difference and p-values for pairwise comparison among classifiers
diffs <- diff(resamps,metric=metric)
diffs.table <- summary(diffs)$table[[1]]

pdf("Plots_Feb2013/Figure_auc_diffs_top.pdf")
levelplot(diffs, what = "differences",region=TRUE,col.regions = rev(heat.colors(20)))
dev.off()

pdf("Plots_Feb2013/Figure_auc_diffs_table_top.pdf")
textplot(diffs.table)
dev.off()




