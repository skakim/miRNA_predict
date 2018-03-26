# Compare classifiers in terms of confidence intervals (95%) over AUC scores
# Mariana Mendoza
# February, 2012

#setwd("/Users/marimendoza/Doutorado/UFRGS/Pesquisas/miRNA/RFMirTarget/")
#source("Scripts/functions.R")
load("Scripts/rf_model_34feat.RData")


library(caret)

# plot accuracy, sensitiviy ans specificity in the same graph, as lines
pdf("Plots_Feb2013/Performance_34feat.pdf")
plot(rf.performance.rfs$values[,1],type="l",ylim=c(70,95),lwd=2,xlab="Features",ylab="Performance")
lines(rf.performance.rfs$values[,2],col="red",lwd=2)
lines(rf.performance.rfs$values[,3],col="blue",lwd=2)
legend("bottomright",legend=c("ACC","SPE","SEN"),col=c("black","red","blue"),inset=0.1,lwd=2)
dev.off()

# plot mcc
pdf("Plots_Feb2013/Performance_34feat2.pdf")
plot(rf.performance.rfs$values[,4],type="l",lwd=2,col="purple",xlab="Features",ylab="Performance")
legend("bottomright",legend=c("MCC"),col=c("purple"),inset=0.1,lwd=2)
dev.off()


