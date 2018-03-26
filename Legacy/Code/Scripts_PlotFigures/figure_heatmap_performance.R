# Compare classifiers in terms of classification probabilities
# Mariana Mendoza
# February, 2012

#setwd("/Users/marimendoza/Dropbox/UFRGS/Pesquisas/miRNA/RFMirTarget/")
source("Scripts/functions.R")
source("Scripts/loadData.R")
load("Scripts/rf_model_34feat.RData")
load("Scripts/rf_model_topfeat.RData")
load("Scripts/rf_model_compare.RData")
load("Scripts/rf_model_compare_topFeat.RData")


library(caret)


# PREDICTION BASED ON ALL FEATURES
rf.predict.test <- predict(rf.fit,data.test.descr,type="prob")
svm.predict.test <- predict(svm.fit,data.test.descr,type="prob")
knn.predict.test <- predict(knn.fit,data.test.descr,type="prob")
j48.predict.test <- predict(j48.fit,data.test.descr,type="prob")
nb.predict.test <- predict(nb.fit,data.test.descr,type="prob")
glm.predict.test <- predict(glm.fit,data.test.descr,type="prob")


# PREDICTION BASED ON TOP FEATURES
rf.top.predict.test <- predict(rf.top.fit,data.test.descr[,colnames(data.train.descr.topFeat)],type="prob")
svm.top.predict.test <- predict(svm.top.fit,data.test.descr[,colnames(data.train.descr.topFeat)],type="prob")
knn.top.predict.test <- predict(knn.top.fit,data.test.descr[,colnames(data.train.descr.topFeat)],type="prob")
j48.top.predict.test <- predict(j48.top.fit,data.test.descr[,colnames(data.train.descr.topFeat)],type="prob")
nb.top.predict.test <- predict(nb.top.fit,data.test.descr[,colnames(data.train.descr.topFeat)],type="prob")
glm.top.predict.test <- predict(glm.top.fit,data.test.descr[,colnames(data.train.descr.topFeat)],type="prob")


all.pred.prob <- cbind(paste("inst",row.names(rf.predict.test),sep=""),rf.predict.test,svm.predict.test,knn.predict.test,j48.predict.test,nb.predict.test,glm.predict.test, data.test.class)

idx.pos <- sample(c(1:length(which(all.pred.prob[,14] == "Target"))),50)
idx.neg <- sample(c(1:length(which(all.pred.prob[,14] == "NonTarget"))),30)


all.pred.true <- all.pred.prob[idx.pos,c(1,which(colnames(all.pred.prob) == "Target"),14)]
colnames(all.pred.true) <- c('ID','RF','SVM','KNN','J48','NB','GLM','CLASS')
all.pred.true.m <- melt(all.pred.true)
# > length(which(all.pred.true[,2]<=0.5))
# [1] 6
# > length(which(all.pred.true[,3]<=0.5))
# [1] 6
# > length(which(all.pred.true[,4]<=0.5))
# [1] 5
# > length(which(all.pred.true[,5]<=0.5))
# [1] 7
# > length(which(all.pred.true[,6]<=0.5))
# [1] 1

all.pred.false <- all.pred.prob[idx.neg,c(1,which(colnames(all.pred.prob) == "NonTarget"),14)]
colnames(all.pred.false) <- c('ID','RF','SVM','KNN','J48','NB','GLM','CLASS')
all.pred.false.m <- melt(all.pred.false)


## heatmap
require(ggplot2)
pdf("Plots_May2013/heatmap_prob_Target2.pdf")
colnames(all.pred.true.m) <- c("obs","class","method","prob")
ggplot(all.pred.true.m, aes(all.pred.true.m$method,all.pred.true.m$obs)) + geom_tile(aes(fill = all.pred.true.m$prob),colour = "black") + scale_fill_gradient(low="white", high="blue")
dev.off()


pdf("Plots_May2013/heatmap_prob_NonTarget2.pdf")
colnames(all.pred.false.m) <- c("obs","class","method","prob")
ggplot(all.pred.false.m, aes(all.pred.false.m$method,all.pred.false.m$obs)) + geom_tile(aes(fill = all.pred.false.m$prob),colour = "black") + scale_fill_gradient(low="white", high="blue")
dev.off()


# ## pairs-plot
# require(GGally)
# ggpairs(all.pred.true, colour='RunOvers', alpha=0.4)
# 
# ## heatmap
# require(ggplot2)
# all.pred.true.m <- melt(cbind(rownames(all.pred.true[1:100,]),all.pred.true[1:100,]))
# colnames(all.pred.true.m) <- c("obs","class","method","prob")
# (q <- ggplot(subset(all.pred.true.m, class=="Yes"), aes(method,obs)) + geom_tile(aes(fill = prob),
#                                                                                  colour = "black") + scale_fill_gradient(low = "white",
#             