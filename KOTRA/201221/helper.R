## helper functions
## X : list, predicted probabilities per each k-fold CV set
## Y : list, true binary responses per each k-fold CV set

#############################################################################################################
## binResids
## binned residual plot for binary responses
## plot predicted probabilities against the standardized residuals against the real value (0, 1) per each bin
## plotted with TRAIN data set, for model diagnostic purpose
## see p.157 of BDA 3ed (Gelman et al.)
## nBins : number of bins to plot, by default sqrt(n)
## ... : additional arguments for plot function
binResids = function(X, Y, nBins, ...){
  old_pars <- par(no.readonly = TRUE)
  on.exit(par(old_pars), add = TRUE)
  
  AVGEXP = AVGRES = vector("list", length = length(X))
  
  ## sort X, Y
  for(i in seq_along(X)){
    ## create bins
    N = length(X[[i]])
    if(missing(nBins)) nBins= ceiling(sqrt(N))
    binsIdx = rep(1:nBins, each = N %/% nBins)
    if(length(binsIdx) < N) binsIdx = c(binsIdx, rep(tail(binsIdx,1), N - length(binsIdx)))
    
    ## set up placeholders
    x = y = sortIdx = xSort = ySort = rep(0, N)
    avgExp = avgRes = rep(0, nBins)
    
    x = X[[i]]
    y = Y[[i]]
    sortIdx = order(x)
    xSort = x[sortIdx]
    ySort = y[sortIdx]
    for(j in 1:nBins){
      avgExp[j] = mean(xSort[binsIdx == j])
      avgRes[j] = mean(ySort[binsIdx == j] - xSort[binsIdx == j]) / sd(ySort[binsIdx ==j] - xSort[binsIdx == j])
    }
    AVGEXP[[i]] = avgExp
    AVGRES[[i]] = avgRes
  }
  
  cutoff = max(3, max(abs(unlist(AVGRES)), na.rm=T))
  if(!hasArg(ylim)) ylim = c(-cutoff, cutoff)
  
  
  par(mfrow=c(1,2))
  plot(avgExp, avgRes, xlim=c(0,1), ...,
       type="n",
       xlab = "Average Predicteds per bins", ylab = "Average Scaled Residuals per bins",
       main = paste0("Binnded Residuals Plot"))
  for(i in seq_along(X)){
    points(AVGEXP[[i]], AVGRES[[i]], pch=1)
  }
  abline(h=1.96, lty=2)
  abline(h=-1.96, lty=2)
  abline(h = 0) 
  
  hist(unlist(AVGRES), nclass= 20,
       xlab = "residuals", prob=T, main= "Histogram of Residuals")
}
#############################################################################################################
# confusion matrix
ConfusionMat = function(predprob,true,cutoff){
  tbl = matrix(0,nrow=2,ncol=2)
  len = length(predprob)
  for(i in 1:len){
    pred = ifelse(predprob[[i]] > cutoff,1,0)
    tbl = tbl + table(pred,true[[i]])
  }
  return(round(tbl/len))
}

#############################################################################################################
## CV ROC and AUC
## https://cran.r-project.org/web/packages/cvAUC/cvAUC.pdf
library(cvAUC)
library(ROCR)

cvROC = function(X, Y){
  old_pars <- par(no.readonly = TRUE)
  on.exit(par(old_pars), add = TRUE)
  par(mfrow=c(1,1))
  out <- cvAUC::cvAUC(X, Y)

  plot(out$perf, col = "grey82", lty = 3, 
                   main = paste0(length(X), "-fold CV ROC"))
  plot(out$perf, col ="red", avg = "vertical", add = TRUE)
  
  tmp = cvAUC::ci.cvAUC(X, Y)
  text(0.6, 0.2, paste0("CV AUC : ", 
                        signif(tmp$cvAUC, 4), 
                        "\n", "CI : ", 
                        signif(tmp$ci[1], 4), " , ", 
                        signif(tmp$ci[2], 4)),
       font =2, adj = 0)
}


cvcutoff = function(pred,true,idx){
  predprob = pred # prob of label 1
  labels = as.factor(true)
  
  predlabel = prediction(predprob,labels)
  perf = performance(predlabel,"acc")
  if(idx ==1)  plot(perf,main="Acc plot",lty=2,col='gray',ylim=c(0,1))
  else plot(perf,main="Acc plot",add=T,lty=2,col='gray')
  
  yloc = max(perf@y.values[[1]])
  xloc = mean(as.numeric(perf@x.values[[1]][perf@y.values[[1]]==yloc]))

  return(xloc)
}

cvACC =function(X, Y){
  par(mfrow=c(1,1))
  xlocs = rep(NA,length(X))
  plot.new()
  for(i in 1:length(X)){
    cutoff = cvcutoff(X[[i]],Y[[i]],idx=i)
    xlocs[i] = cutoff[[1]]
  }
  
  abline(v=mean(xlocs),lty=2,col="red")
  text(mean(xlocs)+0.08,0.05,paste0("Optimal cutoff: ",round(mean(xlocs),2)))
}
#############################################################################################################




#############################################################################################################
## varImp
library(caret)
library(ggplot2)

plotVI = function(tmp){
  # old_pars <- par(no.readonly = TRUE)
  # on.exit(par(old_pars), add = TRUE)
  # par(mar=c(3,8,1,2))
  idx = order(tmp$Overall)
  barplot(tmp$Overall[idx]/ sum(tmp$Overall),
          names = rownames(tmp)[idx],
          horiz= T,
          las = 1,
          col = "skyblue",
          main = "Variable Importance")
}

plotVIlist = function(VIlist, p = 3){
  old_pars <- par(no.readonly = TRUE)
  par(mar=c(3,8,1,2), mfrow=c(1,p))
  on.exit(par(old_pars), add = TRUE)
  for(i in sample(length(VIlist), p)){
    plotVI(VIlist[[i]])
  }
}



#############################################################################################################
upsampling = function(dataset){
  num_up = table(dataset$y)[2] - table(dataset$y)[1]
  samp_idx =1:table(dataset$y)[1]
  samp_idx = sample(samp_idx,num_up,replace=TRUE)
  up_sampled = dataset[dataset$y == 0,][samp_idx,]
  
  table(dataset$y)
  final_case_up = rbind(dataset,up_sampled)
  
  cat("before: ",table(dataset$y),"\n")
  cat("after: ",table(final_case_up$y),"\n")
  return(final_case_up)
}
#############################################################################################################

## fit kcv
library(gbm)
library(randomForest)
FITRES = function(dat, kCV, method = "glm", upsamp=T){
  if(!(method %in% c("glm", "ada", "rf"))) stop("invalid methods")
  
  ## setup placeholders
  FITs = yPRED = yTRUE = vIMP = vector("list", length = length(kCV))
  
  for(i in seq_along(kCV)){
    cat(paste0(i,"th iteration..."),"\n")
    ## split data
    testIdx = kCV[[i]]
    trainIdx = (1:nrow(dat))[-testIdx]
    XtestDf = dat[-trainIdx, colnames(dat) != "y"]
    ytest = dat[-trainIdx, colnames(dat) == "y", drop=T]
    
    datTrain = dat[trainIdx,]
    if(upsamp){ datTrain = upsampling(datTrain) }
    XtrainDf = datTrain[, colnames(datTrain) != "y"]
    ytrain = datTrain[, colnames(datTrain) == "y", drop=T]
    
    if(method == "glm"){
      fit = glm(y ~ ., data = datTrain, family = binomial("logit"))
      vImp = varImp(fit)
    } else if (method == "ada") {
      fit = gbm(y~., distribution = "adaboost", data = datTrain, n.trees = 500)
      vImp = varImp(fit, numTrees = fit$n.trees)
    } else if (method == "rf") {
      fit = randomForest(as.factor(y)~., data = datTrain) # y(factor) : classification is assumed.
      vImp = varImp(fit)
    }
    predprob = predict(fit, newdata = XtestDf, type="prob")
    
    FITs[[i]] = fit
    yPRED[[i]] = predprob[,2] # prob of label 1
    yTRUE[[i]] = ytest
    vIMP[[i]] = vImp
  }
  
  return(list(FITs = FITs, yPRED = yPRED, yTRUE = yTRUE, vIMP = vIMP))
}
#############################################################################################################
# convert factor list to numeric list
conv_lst = function(factor_lst){
  len = length(factor_lst)
  for(i in 1:len){
    factor_lst[[i]] = as.numeric(unlist(factor_lst[i])) - 1
  }
  return(factor_lst)
}


#############################################################################################################



















































