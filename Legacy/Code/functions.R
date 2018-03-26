CompareModelsConfusionMatrix <- function(models,data.descr,data.class) {
  # Computes the confusion matrix and the accuracy, sensititvy and specificity 
  # for several models obtained with train() function from carte package, based on 
  # an independent test set
  #
  # Args:
  #   models: list with models, elements properly named with models' name
  #   data.descr: test data
  #   data.label: class of instances of test data (must have the same number 
  #               of rows than data.descr)
  #
  # Returns:
  #   A list with accuracy, specificity and sensitivity values for all models
  # 
  # Mariana Mendoza, February 2012
  #
  stopifnot(dim(data.descr)[1] == length(data.class))
  
  if(missing(models)) {
    stop("A list with the models (return by train() function from caret package) must be given as parameter for this function.")
  }

  spe <- NULL
  sen <- NULL
  acc <- NULL  
  methods.names <- names(models)

  for (ii in 1:length(methods.names)) {
    predict.test <- predict(models[[methods.names[ii]]],data.descr,type="raw")
    perf <- confusionMatrix(predict.test,data.class,positive="Target")
    spe <- c(spe, perf$byClass["Specificity"])
    sen <- c(sen, perf$byClass["Sensitivity"])
    acc <- c(acc, perf$overall["Accuracy"])
  }
  names(sen) <- methods.names
  names(spe) <- methods.names
  names(acc) <- methods.names
  

  results <- list(ACC=acc,SPE=spe,SEN=sen)
  return(results)
}

#---

CompareModelsAUCScore <- function(models,data.descr,data.class,plotROC=FALSE,filename="roc_curve.pdf") {
  # Computes AUC scores for several models obtained with train() function from caret 
  # package based on an independent test set
  #
  # Args:
  #   models: list with models, elements properly named with models' name
  #   data.descr: test data
  #   data.label: class of instances of test data (must have the same number 
  #               of rows than data.descr)
  #   plotROC: must ROC curve be ploted? (default:FALSE) 
  #   filename: in case of plotROC=TRUE, name of the output PDF file
  #
  # Returns:
  #   A vector with the AUC scores of all models
  # 
  # Mariana Mendoza, February 2012
  #
  stopifnot(dim(data.descr)[1] == length(data.class))
  
  if(missing(models)) {
    stop("A list with the models (return by train() function from caret package) must be given as parameter for this function.")
  }
  
  auc <- NULL
  
  methods.names <- names(models)
  color <- rainbow(length(methods.names))   
  
  if (plotROC) { pdf(filename)}
  for (ii in 1:length(methods.names)) {
    predict.test <- predict(models[[methods.names[ii]]],data.descr,type="prob")
    predict.test.pred <- prediction(predict.test[,"Target"],data.class)
    perf <- performance(predict.test.pred,"auc")
    auc <- c(auc,perf@y.values[[1]])
    if (plotROC) {
      perf <- performance(predict.test.pred,"tpr","fpr")
      if (ii == 1) {
        plot(perf,col=color[ii],lwd=2) 
      } else {
        plot(perf,col=color[ii],lwd=2,add=TRUE)  
      }
    }
  }
  names(auc) <- methods.names
  if (plotROC) {
    abline(0,1,col="black",lwd=1)
    legend("bottomright",legend=methods.names,lwd=2,col=color)
  }
  if (plotROC) { dev.off()}
  return(auc)
}

#---

ComputeFeaturesImportance <- function(model,criteria) {
  # Computes estimative of feature importance based on the OOB data
  #
  # Args:
  #   model: model trained with the randomForest() function
  #   criteria: criteria for feature importance, "accuracy" or "gini". Default is "accuracy". 
  #
  # Returns:
  #   A vector with the accuracy, specificity, sensitivity and Matthews Correlation Coefficient computed based on 
  # the OOB data.
  # 
  # Mariana Mendoza, November 2012
  #
  if(missing(model)) {
    stop("A randomForest class object must be given as parameter for this function.")
  }
  
  if(missing(criteria)) {
    feat.importance <- sort((model$importance[,"MeanDecreaseAccuracy"])*100,decreasing = TRUE)
  }
  else{
    if((criteria != "accuracy") && (criteria != "gini")) {
      stop("Criteria for features importance must be either 'accuracy' or 'gini'.")
    }
    else {
      if(criteria == "accuracy") {
        feat.importance <- sort((model$importance[,"MeanDecreaseAccuracy"])*100,decreasing = TRUE)
      }
      else {
        feat.importance <- sort((model$importance[,"MeanDecreaseGini"]),decreasing = TRUE)
      }
    }
  }
  
  feat.importance.df <- as.data.frame(feat.importance)
  
  return(feat.importance.df)
}

#---

ComputeOOBStatistics <- function(model,verbose=FALSE) {
  # Computes the statistics from the OOB data of a random forest model
  #
  # Args:
  #   model: model trained with the randomForest() function
  #
  # Returns:
  #   A vector with the accuracy, specificity, sensitivity and Matthews Correlation Coefficient computed based on 
  # the OOB data.
  # 
  # Mariana Mendoza, November 2012
  #
  if(missing(model)) {
    stop("A randomForest class object must be given as parameter for this function.")
  }
  
  ACC <- ((model$confusion[1,1] + model$confusion[2,2]) / nrow(model$votes))*100
  SPE <- (model$confusion[1,1] / (model$confusion[1,1] + model$confusion[1,2]))*100
  SEN <- (model$confusion[2,2] / (model$confusion[2,1] + model$confusion[2,2]))*100
  MCC <- MatthewsCorrelationCoefficient(model$confusion)
  
  performance <- cbind(ACC,SPE,SEN,MCC)
  colnames(performance) <- c("ACC","SPE","SEN","MCC")
  
  if (verbose) {
    print(performance)
  }
  
  return(performance)
}

#---

Confusion2Performance <- function(confusionMat,verbose=FALSE) {
  # Extract model performance from confusionMatrix structure returned by caret
  # package function confusionMatrix(). 
  #
  # MCC, specificity, sensitivity, accuracy and confidence levels for accuracy
  # are returned in a list. 
  #
  # Args:
  #   confusionMat: confusion matrix returned by confusionMatrix() function from
  # caret package
  #
  # Returns:
  #   A list with elements 'prediction' and 'accuracy'
  # 
  # Mariana Mendoza, November 2012
  #
  if(missing(confusionMat)) {
    stop("A confusioMatrix class object must be given as parameter for this function.")
  }
  
  if (verbose) {
    print(confusionMat$table)
  }
  
  MCC <- MatthewsCorrelationCoefficient(confusionMat$table)
  SEN <- confusionMat$byClass[1]*100
  SPE <- confusionMat$byClass[2]*100
  ACC <- confusionMat$overall[1]*100
  ACC.low <- confusionMat$overall[3]*100
  ACC.up <- confusionMat$overall[4]*100

  performance <- cbind(ACC,SPE,SEN,MCC)
  colnames(performance) <- c("ACC","SPE","SEN","MCC")
  row.names(performance) <- NULL
  
  accuracy.ci <- cbind (ACC.low, ACC.up)
  colnames(accuracy.ci) <- c("AccuracyLower","AccuracyUpper")
  row.names(accuracy.ci) <- NULL

  
  if (verbose) {
    print(performance)
    print(accuracy.ci)
  }
  
  output <- list(performance,accuracy.ci)
  names(output)<-c("performance","accuracy")

  return(output)
}

#---


CreatePredictionsListCV <- function(oof,repeats) {
  # Creates a list with the predicted and observed values for each CV repetition
  # Args:
  #   oof: predictions for the out-of-fold data (.$pred element from the object 
  # returned by train() function when option savePredictions=TRUE)
  #   repeats: number of repeats of the cross-validation
  #
  # Returns:
  #   A list with the predictions and the labels, with size equal the number of CV repeatitions
  # 
  # Mariana Mendoza,February 2012
  #
  if(missing(oof)) {
    stop("The predictions for the out-of-fold data must be given as parameter for this function.")
  }
  
  if(missing(repeats)) {
    stop("It is necessary to inform the number of repetitions for the cross-validation algorithm")
  }
  
  folds <- unique(oof$Resample)
  predictions <- NULL
  labels <- NULL
  
  #for each repeat, build vector of predicted and observed values
  for(ii in 1:repeats) {
    resample <- oof[is.element(oof$Resample,folds[grep(paste(".Rep",ii,sep=""),folds)]),]
    
    predictions <- c(predictions,list(as.matrix(resample["Target"])))
    labels <- c(labels,list(as.matrix(resample["obs"])))
  }
  
  results <- list(predictions,labels)
  names(results) <- c("Predictions","Labels")
  return(results)
}

#---

EstimatePerformanceCV <- function(oof,repeats) {
  # Computes the average confusion matrix and the standard deviation from the output
  # of a (repeated) cross-validation
  #
  # Args:
  #   oof: predictions for the out-of-fold data (.$pred element from the object 
  # returned by train() function when option savePredictions=TRUE)
  #   repeats: number of repeats of the cross-validation
  #
  # Returns:
  #   A list with the average confusion matrix and the standard deviation when repeated
  # cross-validation is performed.
  # 
  # Mariana Mendoza, December 2012
  #
  if(missing(oof)) {
    stop("The predictions for the out-of-fold data must be given as parameter for this function.")
  }
  
  if(missing(repeats)) {
    stop("It is necessary to inform the number of repetitions for the cross-validation algorithm")
  }
  
  folds <- unique(oof$Resample)
  all.confusion <- NULL
  all.performance <- NULL
  
  for(ii in 1:repeats) {
    resamp <- oof[is.element(oof$Resample,folds[grep(paste(".Rep",ii,sep=""),folds)]),]
    oof.conf <- confusionMatrix(resamp$pred,resamp$obs,positive="Target")
    all.confusion <- c(all.confusion,list(oof.conf$table))
    oof.performance <- Confusion2Performance(oof.conf)
    #print(oof.conf)
    #print(oof.performance)
    all.performance <- c(all.performance,list(oof.performance$performance))

  }
  names(all.confusion) <- paste("Rep",c(1:repeats),sep="")
  
  TP <- NULL
  FP <- NULL
  TN <- NULL
  FN <- NULL
  
  ACC <- NULL
  SPE <- NULL
  SEN <- NULL
  MCC <- NULL
  
  for (ii in 1:length(all.confusion)) {
    TP <- c(TP, all.confusion[[ii]]["Target","Target"])
    FN <- c(FN, all.confusion[[ii]]["NonTarget","Target"])
    TN <- c(TN, all.confusion[[ii]]["NonTarget","NonTarget"])
    FP <- c(FP, all.confusion[[ii]]["Target","NonTarget"])
    
    ACC <- c(ACC,all.performance[[ii]][,"ACC"])
    SEN <- c(SEN,all.performance[[ii]][,"SEN"])
    SPE <- c(SPE,all.performance[[ii]][,"SPE"])
    MCC <- c(MCC,all.performance[[ii]][,"MCC"])
  }
  
  # compute mean and deviation for confusion matrices
  TP.mean <- mean(TP)
  TP.sd <- sd(TP)
  if (is.na(TP.sd)) {TP.sd <- 0}
  
  FP.mean <- mean(FP)
  FP.sd <- sd(FP)
  if (is.na(FP.sd)) {FP.sd <- 0}
  
  TN.mean <- mean(TN)
  TN.sd <- sd(TN)
  if (is.na(TN.sd)) {TN.sd <- 0}
  
  FN.mean <- mean(FN)
  FN.sd <- sd(FN)
  if (is.na(FN.sd)) {FN.sd <- 0}

  mean.confusion <- matrix(c(TN.mean,FN.mean,FP.mean,TP.mean),byrow=TRUE,nrow=2,ncol=2,dimnames=list(c("Pred NonTarget","Pred Target"),c("Ref NonTarget","Ref Target")))
  sd.confusion <- matrix(c(TN.sd,FN.sd,FP.sd,TP.sd),byrow=TRUE,nrow=2,ncol=2,dimnames=list(c("Pred NonTarget","Pred Target"),c("Ref NonTarget","Ref Target"))) 
  
  # compute mean and deviation for performance metrics
  ACC.mean <- mean(ACC)
  ACC.sd <- sd(ACC)
  SPE.mean <- mean(SPE)
  SPE.sd <- sd(SPE)
  SEN.mean <- mean(SEN)
  SEN.sd <- sd(SEN)
  MCC.mean <- mean(MCC)
  MCC.sd <- sd(MCC)
  
  mean.metrics <- c(ACC.mean, SPE.mean, SEN.mean, MCC.mean)
  names(mean.metrics) <- c("ACC","SPE","SEN","MCC")
  sd.metrics <- c(ACC.sd, SPE.sd, SEN.sd, MCC.sd)
  names(sd.metrics) <- c("ACC","SPE","SEN","MCC")
  
  performance.cv <- list(mean.confusion,sd.confusion,mean.metrics,sd.metrics)
  names(performance.cv) <- c("mean.confusion","deviation.confusion","mean.metrics","deviation.metrics")
  return(performance.cv)
}

#---


EstimateErrorRateCV <- function(oof,repeats) {
  # Computes the average error rate for each class based on (repeated) cross-validation
  #
  # Args:
  #   oof: predictions for the out-of-fold data (.$pred element from the object 
  # returned by train() function when option savePredictions=TRUE)
  #   repeats: number of repeats of the cross-validation
  #
  # Returns:
  #   A list with the average confusion matrix and the standard deviation when repeated
  # cross-validation is performed.
  # 
  # Mariana Mendoza, February 2012
  #
  if(missing(oof)) {
    stop("The predictions for the out-of-fold data must be given as parameter for this function.")
  }
  
  if(missing(repeats)) {
    stop("It is necessary to inform the number of repetitions for the cross-validation algorithm")
  }
  
  folds <- unique(oof$Resample)
  all.confusion <- NULL
  
  
  for(ii in 1:repeats) {
    resamp <- oof[is.element(oof$Resample,folds[grep(paste(".Rep",ii,sep=""),folds)]),]
    oof.conf <- confusionMatrix(resamp$pred,resamp$obs,positive="Target")
    all.confusion <- c(all.confusion,list(oof.conf$table))
    #oof.performance <- Confusion2Performance(oof.conf)
    #print(oof.conf)
    #print(oof.performance)
    #all.performance <- c(all.performance,list(oof.performance$performance))
    
  }
  names(all.confusion) <- paste("Rep",c(1:repeats),sep="")
  
  pos.error <- NULL
  neg.error <- NULL
  
  for (ii in 1:length(all.confusion)) {
    totalExamples <- colSums(all.confusion[[ii]])
    names(totalExamples) <- colnames(all.confusion[[ii]])
    FN <- all.confusion[[ii]]["NonTarget","Target"]
    FP <- all.confusion[[ii]]["Target","NonTarget"]
    pos.error <- c(pos.error,(FN/totalExamples["Target"])*100)
    neg.error <- c(neg.error,(FP/totalExamples["NonTarget"])*100)
  }
  
  
  error.cv <- list(mean(pos.error),sd(pos.error),mean(neg.error),sd(neg.error))
  names(error.cv) <- c("mean.positive","deviation.positive","mean.negative","deviation.negative")
  return(error.cv)
}

#---

MatthewsCorrelationCoefficient <- function(confusionMat) {
  # Compute Matthews Correlation Coefficient
  #
  # Args:
  #   confMat: confusion matrix
  #
  # Returns:
  #   MCC 
  #
  TN <- confusionMat[1,1]
  FN <- confusionMat[1,2]
  FP <- confusionMat[2,1]
  TP <- confusionMat[2,2]
  
  num <- as.numeric((TP * TN) - (FP * FN))
  den1 <- as.numeric((TP + FP) * (TP + FN))
  den2 <- as.numeric((TN + FP) * (TN + FN))
  
  MCC <- num/sqrt(den1*den2)
  
}

#---

PermutationTest <- function(model,data.class,auc,N=1000) {
  # Permute class labels N times and evaluate prediction. 
  # Compute p-value, which represents the fraction of randomized samples where the
  # classifier behaved better in the random data than in the original data. 
  # Intuitively, it measures how likely the observed accuracy would be obtained by chance.  
  #
  # Args:
  #   model: prediction, obtained by predict() function with type="prob" option
  #   data.class: true class of data
  #   auc: auc score obtained form evaluation on true classe values
  #   N: number of permutaiton (default = 1000)
  #
  # Returns:
  #   p-value
  #
  options(digits=10)
  
  is.higher <- 0
  for (ii in 1:N) {
    permutation <- permute(c(1:length(data.class)))
    o <- order(permutation)
    permutation.class <- data.class[o]
    pred <- prediction(model[,"Target"],permutation.class)
    perf <- performance(pred,"auc")
    if (perf@y.values[[1]] > as.numeric(auc)) {
      is.higher <- is.higher + 1
    }
  }
  pvalue <- is.higher/N
  return(pvalue)
}

# ---

PlotDensityCurves <- function(methods,metric,resamples,rug=FALSE) {
  # Computes density curves and rug for all methods in a single graph
  #
  # Args:
  #   model: character vetor with methods names (following colnames of resamples$values variable)
  #   metric: name of metric to use: "ROC", "Spec", "Sens"
  #   resamples: 'values' variable obtained with caret resamples() function (resamps$values)
  #
  # 
  # Mariana Mendoza, February 2012
  #
  numMethods <- length(methods)
  for (ii in 1:numMethods) {
    if (length(grep(methods[[ii]],colnames(resamples))) == 0) {
      stop("Model names informed must match the same methods used by resample() function.")
    }
  }  
  
  index <- paste(methods,"~",metric,sep="")
  
  #compute x limits
  xlimit.min<-round(min(resamples[index]),digits=1)-0.1
  xlimit.max<-round(max(resamples[index]),digits=1)+0.1
  
  #compute y limits
  ylimit.max <- 0
  for(ii in 1:numMethods) {
    d<-density(as.numeric(unlist(resamples[index[ii]])))
    if(max(d$y) > ylimit.max) {ylimit.max <- max(d$y)}
  }
  ylimit.max <- ceiling(ylimit.max)
  
  color <- rainbow(numMethods)   
  
  plot(x=NULL,y=NULL,xlim=c(xlimit.min,xlimit.max),ylim=c(0,ylimit.max),xlab="AUC", ylab="Density")
  for(ii in 1:numMethods) {
    lines(density(as.numeric(unlist(resamples[index[ii]]))),lwd=2,col=color[ii]) 
    if(rug) { 
      rug(as.numeric(unlist(resamples[index[ii]])),lwd=1.5,col=color[ii]) 
    } else {
      points(as.numeric(unlist(resamples[index[ii]])),y=rep(0,dim(resamples[index[ii]])[1]),lwd=1,col=color[ii])
    }
    
  }
  legend("topleft",legend=methods,col=color,lwd=1)
}

#---

RestrictedForwardSelection <- function(data,model,importance,verbose=FALSE) {
  # Performs feature selection by training several RF models adding one feature at a time,
  # following a decreasing order of importance.
  #
  # Args:
  #   data: training data
  #   model: model trained with the randomForest() function
  #   importance: estimative of features importance returned by ComputeFeaturesImportance()
  #   verbose: print the features list if TRUE; don't print if FALSE
  #
  # 
  # Returns:
  #   A matrix (dimension #features x 4) with the accuracy, specificity, sensitivity and Matthews 
  #   Correlation Coefficient computed based on the OOB data for each iteration.
  # 
  # Mariana Mendoza, November 2012
  #
  features <- NULL
  performance.RFS <- NULL
  num.features <- length(colnames(data))
  
  for (ii in 1:(num.features-1)) {
    features <- c(features,which(colnames(data) == row.names(importance)[ii]))
    
    if (verbose) print(features)
    
    model <- TrainModel(data[,c(features,num.features)],seed=2)
    performance <- ComputeOOBStatistics(model,verbose)
    performance.RFS <- rbind(performance.RFS,performance)
  }
  
  result <- list(performance.RFS,features)
  names(result) <- c("values","features")
  return(result)
}

#---

TestModel <- function(model,data,verbose=TRUE) {
  # Test a trained random forest model using and independent test set
  #
  # Args:
  #   model: model trained with the randomForest() function
  #   data: independent test set
  #
  # Returns:
  #   A vector with the accuracy, specificity, sensitivity and Matthews 
  #   Correlation Coefficient computed based on the given independent data.
  #
  # Mariana Mendoza, November 2012
  #
  if(missing(data)) {
    stop("User needs to indicate a data frame with testing data.")
  }
  
  if(!is.element("Class", colnames(data))){ 
    stop("A 'Class' column with the correct label of each instance is mandatory.")
  }
  
  library(caret)
  
  pred <- predict(model, data)
  pred.conf <- confusionMatrix(data$Class,pred,positive="Target")
  output <- Confusion2Performance(pred.conf,verbose=verbose)
  
  output <- c(list(pred),output)
  names(output)<-c("prediction","performance","accuracy")
  
  return(output)
}

#---

TrainRFModel <- function(data,numTrees,numVar,seed,weights) {
  # Train a random forest model based on the randomForest R package
  #
  # Args:
  #   data: training data
  #   numTrees: number of trees to grow (optional)
  #   numVar: number of predictors (optional)
  #   seed: seed to repeat experiments (optional)
  #   weights: class weights for missclassification 
  #
  # Returns:
  #   An object of class randomForest
  #
  # Mariana Mendoza, November 2012
  #
  if(!is.element("Class", colnames(data))){ 
    stop("A 'Class' column with the correct label of each instance is mandatory.")
  }
  
  if(missing(numTrees)) {
    numTrees <- 500
  }
  
  if(missing(numVar)) {
    numVar <- floor(sqrt(ncol(data)-1))
  }
  
  if(missing(seed)) {
    seed <- 2
  }
  
  if(missing(weights)) {
    weights <- NULL
  }
  
  library(randomForest)
  
  set.seed(seed)
  
  model <- randomForest(Class ~ ., data, importance=TRUE, ntree=numTrees, mtry=numVar,weights=weights)

  return(model)
}

#---

TuneModel <- function(data,labels,method,metric,myFolds,number,repeats,myGrid,seed,weights,verbose=TRUE) {
  # Tune the 'mtry' parameters and train a random forest model based on the caret
  # R package, using repeated cross-validation to estimate the model performance.
  #
  # Args:
  #   data: training data
  #   labels: correct reponse for each instance in training data
  #   method: method to be run, based on 'caret' package description
  #   metric: metric to maximize. "Accuracy" (default) or "Kappa"
  #   myFolds: folds configuration created by createMultiFolds() function (optional)
  #   number: number of folds for cross-validation
  #   repeats: number of times to repeat X-fold cross-validation
  #   myGrid: used in parameters optimization, created by expand.grid() function
  #   seed: seed to repeat experiments (optional)
  #   weights: class weights for missclassification (optional)
  
  #
  # Returns:
  #   A list, as specified by the 'caret' package
  #
  # Mariana Mendoza, November 2012
  #  
  if(missing(data)) {
    stop("User needs to indicate a data frame with training data.")
  }
  
  if(missing(metric)) {
    metric = "Accuracy"  
  }
  
  if(missing(seed)) {
    seed <- 9
  }
  
  if(missing(weights)) {
    weights<-NULL
  }
  
  set.seed(seed)
  
  library('caret')

  if(missing(myFolds)) {
    myControl <- trainControl(## 10-fold CV
        method = "repeatedcv",
        ## fold
        number = number,
        ## repetitions
        repeats = repeats,
        savePredictions=TRUE,
        classProbs = TRUE,
        summaryFunction = twoClassSummary)
  } else {
    myControl <- trainControl(## 10-fold CV
        method = "repeatedcv",
        number = number,
        repeats = repeats,
        index = myFolds,
        savePredictions=TRUE,        
        classProbs = TRUE,
        summaryFunction = twoClassSummary)
  }
  
  modelFit <- train(data, labels,
               method = method,
               metric = metric,
               trControl = myControl,
               tuneGrid = myGrid,
               weights = weights,
               selectionFunction = "oneSE")
  
  if (verbose) {
    print(modelFit)
  }
  
  return(modelFit)
}