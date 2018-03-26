setwd("/Users/marimendoza/Dropbox/UFRGS/Pesquisas/miRNA/RFMirTarget/")

library(randomForest)
library(caret)

source("Scripts/functions.R")


#-------------------------------------------------------------------------------
#LOAD DATA
source("loadData.R")

#-------------------------------------------------------------------------------
#TUNE AND TRAIN RF MODEL WITH 34 FEATURES

source("Scripts/main_rf_34feat.R")

#-------------------------------------------------------------------------------
#TUNE AND TRAIN RF MODEL WITH FOR EACH FEATURE CATEGORY

source("Scripts/main_rf_featcats.R")

#-------------------------------------------------------------------------------
#TUNE AND TRAIN RF MODEL WITH TOP FEATURES

source("Scripts/main_rf_topfeat.R")

#-------------------------------------------------------------------------------
#TUNE AND TRAIN OTHER CLASSIFIERS

source("Scripts/main_compare_classifiers.R")
