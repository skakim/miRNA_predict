# Compare classifiers in terms of confidence intervals (95%) over AUC scores
# Mariana Mendoza
# February, 2012

#setwd("/Users/marimendoza/Doutorado/UFRGS/Pesquisas/miRNA/RFMirTarget/")
#source("Scripts/functions.R")
load("Scripts/rf_model_34feat.RData")
load("Scripts/rf_model_topfeat.RData")

library(caret)
library(randomForest)



x<-margin(rf.model,data.test.class)
y<-margin(rf.top.model,data.test.class)
w<-cbind(x,y)
plot(w[,1],w[,2])

targets <- which(row.names(w) == "Target")
nontargets <- which(row.names(w) == "NonTarget")

colors <- c(1:length(w[,1]))
colors[targets] <- "red"
colors[nontargets] <- "blue"



type <- rep(16,length(w[,1]))
type[which(w[,1] <=0)] <- 1
length(which(w[,1] <=0))/length(w[,1])
# [1] 0.1262136

pdf("margin_34featmodel.pdf")
plot(w[,1],col=colors,pch=type,lwd=1.2,ylab = "Margin", xlab = "Instance")
abline(a=0,b=0,lty=2)  
dev.off()

type <- rep(16,length(w[,2]))
type[which(w[,2] <=0)] <- 1
length(which(w[,2] <=0))/length(w[,2])
# [1] 0.1043689

pdf("margin_12featmodel.pdf")
plot(w[,2],col=colors,pch=type,lwd=1.2,ylab = "Margin", xlab = "Instance")
abline(a=0,b=0,lty=2)  
dev.off()

pdf("margin_correlation.pdf")
plot(w[,1],w[,2],col=c(rep("red",length(which(row.names(w) == "Target"))),rep("blue",length(which(row.names(w) == "NonTarget")))),lwd=1.5,xlab = "34 features model",ylab = "top 12 features model")
dev.off()



cor.test(w[targets,1],w[targets,2])
# Pearson's product-moment correlation
# 
# data:  w[targets, 1] and w[targets, 2] 
# t = 52.4452, df = 480, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0 
# 95 percent confidence interval:
#  0.9082595 0.9349824 
# sample estimates:
#       cor 
# 0.9227221 

cor.test(w[nontargets,1],w[nontargets,2])
# Pearson's product-moment correlation
# 
# data:  w[nontargets, 1] and w[nontargets, 2] 
# t = 36.0769, df = 340, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0 
# 95 percent confidence interval:
#  0.8661811 0.9105071 
# sample estimates:
#      cor 
# 0.890437

cor.test(w[,1],w[,2])
# Pearson's product-moment correlation
# 
# data:  w[, 1] and w[, 2] 
# t = 63.5842, df = 822, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0 
# 95 percent confidence interval:
#  0.8993061 0.9224752 
# sample estimates:
#       cor 
# 0.9116119 
