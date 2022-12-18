# rm(list=ls())
# install.packages("party")
# install.packages("ROCR")
setwd("C:/Users/user/Seungjun/KOTRA")
load("case_data.Rdata")
load("case_data_upsampled.Rdata")
library(randomForest)
library(party)
library(ROCR)
# 
# RF_case = randomForest(y~.,data=final_case1,importance=TRUE)
# RF_case$importance
# varImpPlot(RF_case)
# 
# predictions = RF_case$votes[,2] # prob of label 1
# labels = as.factor(final_case1$y)
# 
# 
# pred_val = prediction(predictions,labels)
# pred = performance(pred_val,"tpr","fpr")
# plot(pred,main="ROC plot")
# pred_AUC = performance(pred_val,"auc")
# text(0.5,0.5,paste("AUC = ",format(pred_AUC@y.values[[1]], digits=5, scientific=FALSE)))
# 
# opt.cut = function(perf, pred){
#   cut.ind = mapply(FUN=function(x, y, p){
#     d = (x - 0)^2 + (y-1)^2
#     ind = which(d == min(d))
#     c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
#       cutoff = p[[ind]])
#   }, perf@x.values, perf@y.values, pred@cutoffs)
# }
# res = opt.cut(pred,pred_val)
# print(opt.cut(pred,pred_val))
# # [,1]
# # sensitivity 0.8494624
# # specificity 0.8504673
# # cutoff      0.5014893
# points(1-res[2,1],res[1,1],pch=20,col="red")
# 
# text(1-res[2,1]+0.05,res[1,1],round(res[3,1],4))



#------------------------------------------------------#
# functions
#------------------------------------------------------#
# cutoff 
opt.cut = function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p){
    d = (x - 0)^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}

show_RF = function(data,seed,prop){
  set.seed(seed)
  sn = sample(1:nrow(data),size=nrow(data)*prop)
  case_train = data[sn,]
  case_test = data[-sn,]
  
  # ctree_case = ctree(y~.,data=case_train)
  # plot(ctree_case,type="simple")
  
  RF_case = randomForest(y~.,data=case_train,importance=TRUE) # ,ntree=500,mtry,replace=TRUE
  
  #plotting
  plot(RF_case,main="OOB error")
  varImpPlot(RF_case,main="Feature Importance")
  
  # AUC
  predictions = predict(RF_case,newdata=case_test,type="prob") # prob of label 1
  labels = as.factor(case_test$y)
  pred_val = prediction(predictions[,2],labels)
  pred = performance(pred_val,"tpr","fpr")
  plot(pred,main="ROC plot")
  pred_AUC = performance(pred_val,"auc")
  text(0.5,0.5,paste("AUC = ",format(pred_AUC@y.values[[1]], digits=5, scientific=FALSE)))

  #cutoff
  res = opt.cut(pred,pred_val)
  points(1-res[2,1],res[1,1],pch=20,col="red")
  text(1-res[2,1]+0.1,res[1,1],paste("Optimal cutoff:",round(res[3,1],4)))
  
  # summary 
  cat("Result of Random Forest\n")
  print(RF_case)
  cat("------------------------\n")
  cat("Feature Importance\n")
  print(RF_case$importance)
  cat("------------------------\n")
  # eval
  RF_casep = predict(RF_case,newdata=case_test)
  cat("Predicted\n")
  print(table(case_test$y,RF_casep))
  cat("------------------------\n")
  #test
  train_acc = sum(RF_casep == case_test$y) / nrow(case_test) * 100
  #train
  test_acc = sum(RF_case$predicted == case_train$y) / nrow(case_train) * 100
  cat("train_acc: ",train_acc," test_acc: ",test_acc,"\n")
}


#------------------------------------------------------#
# case1 no upsampling
#------------------------------------------------------#
show_RF(data=final_case1,seed=123,prop=0.7)

# test_acc:  71.85471
# AUC : 0.73083
#------------------------------------------------------#
# case1 upsampling 
#------------------------------------------------------#
show_RF(data=final_case1_up,seed=123,prop=0.7)
# acc는 별로 차이가 없지만 AUC로 봤을 때는 차이가 생긴다.
# test_acc:  70.69267
# AUC : 0.76327
#------------------------------------------------------#
# case2 no upsampling
#------------------------------------------------------#
show_RF(data=final_case2_imp,seed=123,prop=0.7)
# test_acc:  75.60615
# AUC : 0.76774
#------------------------------------------------------#
# case2 upsampling
#------------------------------------------------------#
show_RF(data=final_case2_up,seed=123,prop=0.7)
# acc도 그렇고 auc까지 엄청나게 올라간다..
# test_acc:  84.98088
# AUC : 0.94211
#------------------------------------------------------#
# case4 no upsampling
#------------------------------------------------------#
show_RF(data=final_case4_imp,seed=123,prop=0.7)
# test_acc:  76.45286
# AUC : 0.799
#------------------------------------------------------#
# case4 upsampling
#------------------------------------------------------#
show_RF(data=final_case4_up,seed=123,prop=0.7)
# acc도 그렇고 auc까지 엄청나게 올라간다..
# test_acc:  86.95993
# AUC : 0.96351




