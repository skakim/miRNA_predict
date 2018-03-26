# Compare classifiers in terms of ROC curves
# Mariana Mendoza
# February, 2012

#setwd("/Users/marimendoza/Doutorado/UFRGS/Pesquisas/miRNA/RFMirTarget/")
source("Scripts/functions.R")
load("Scripts/rf_model_34feat.RData")
load("Scripts/rf_model_featcat.RData")
load("Scripts/rf_model_topfeat.RData")

library(ROCR)

numMethods <- 4
cols <- rev(rainbow(numMethods))
ptype <- c(0:numMethods-1)
ltype <- c(2,3,4,5)

#setEPS()
#postscript("ROC_average_rfvscf.eps",height=4,width=6.5)
#o<-par(xpd=T, mar=c(4,4,1.3,1)+0.1) #

pdf("Plots_Feb2013/Figure_roc_average_rf.pdf")

# RF 34 FEATURES
oof <- rf.fit$pred
oof <- oof[oof$.mtry==rf.fit$bestTune[,'.mtry'],]
repeats<-rf.fit$control$repeats
predictions.cv.rf <- CreatePredictionsListCV(oof,repeats)
predictions.rf <- prediction(predictions.cv.rf$Predictions,predictions.cv.rf$Labels)
perf <- performance(predictions.rf,"tpr","fpr")
plot(perf,lwd=1,avg="threshold",type="l",col=cols[1],lty=ltype[1])#,pch=ptype[3])


# RF TOP 12
oof <- rf.top.fit$pred
oof <- oof[oof$.mtry==rf.top.fit$bestTune[,'.mtry'],]
repeats<-rf.top.fit$control$repeats
predictions.cv.rf.top <- CreatePredictionsListCV(oof,repeats)
predictions.rf.top <- prediction(predictions.cv.rf.top$Predictions,predictions.cv.rf.top$Labels)
perf <- performance(predictions.rf.top,"tpr","fpr")
plot(perf,lwd=1,avg="threshold",type="l",col=cols[2],lty=ltype[2],add=TRUE)#,pch=ptype[3])


# RF TOP 10
oof <- rf.nocorr.fit$pred
oof <- oof[oof$.mtry==rf.nocorr.fit$bestTune[,'.mtry'],]
repeats<-rf.nocorr.fit$control$repeats
predictions.cv.rf.nocorr <- CreatePredictionsListCV(oof,repeats)
predictions.rf.nocorr <- prediction(predictions.cv.rf.nocorr$Predictions,predictions.cv.rf.nocorr$Labels)
perf <- performance(predictions.rf.nocorr,"tpr","fpr")
plot(perf,lwd=1,avg="threshold",type="l",col=cols[3],lty=ltype[3],add=TRUE)#,pch=ptype[2],)

# RF TOP 10 weighted
oof <- rf.nocorr.wgt.fit$pred
oof <- oof[oof$.mtry==rf.nocorr.wgt.fit$bestTune[,'.mtry'],]
repeats<-rf.nocorr.wgt.fit$control$repeats
predictions.cv.rf.nocorr.wgt <- CreatePredictionsListCV(oof,repeats)
predictions.rf.nocorr.wgt <- prediction(predictions.cv.rf.nocorr.wgt$Predictions,predictions.cv.rf.nocorr.wgt$Labels)
perf <- performance(predictions.rf.nocorr.wgt,"tpr","fpr")
plot(perf,lwd=1,avg="threshold",type="l",col=cols[4],lty=ltype[4],add=TRUE)#,pch=ptype[3],add=TRUE)



abline(0,1,col="black",lwd=1)

# Draw a legend.
legend(0.6, 0.6, cex=.8, fill=NULL, border="white", col=cols, c('RF all', 'RF Top12', 'RF Top10','RF Top10 weighted'),lwd=2, lty=ltype)

dev.off()

rm(oof,repeats,perf)

