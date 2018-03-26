# Plot error rates for the RF model
# Mariana Mendoza
# February, 2012

setwd("/Users/marimendoza/Doutorado/UFRGS/Pesquisas/miRNA/RFMirTarget/")
source("Scripts/functions.R")
load("Scripts/main_rf_34feat.R")


# plot error rates

# create Line Chart
# get the range for the x and y axis 
xrange <- c(0,500)
yrange <- c(0.1,0.35)

pdf("Plots_Feb2013/errorModelRF2.pdf")
# set up the plot 
plot(xrange, yrange, type="n", xlab="Trees",ylab="Error") 
colors <- rainbow(3) 
linetype <- c(1:3) 
    
# add lines 
for (i in 1:3) { 
  lines(rf.fit$finalModel$err.rate[,i], type="l", lwd=4,lty=linetype[i], col=colors[i])
} 
    
# add a legend 
legend("topright",c("OOB", "Non-Target","Target"), cex=0.8, lwd=3,col=colors, lty=linetype,inset=0.01)
   
# close pdf()
dev.off()
