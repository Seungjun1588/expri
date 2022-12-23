rm(list=ls())
setwd("C:/Users/user/Seungjun/KOTRA/201221")

source("helper.r")
library(tidyverse)
library(openxlsx)
library("xlsx")
#############################################################################################################
## create k-cv folds
dattype = "imp"
upsamp=F
datdir = paste0("data/", "case_data_", dattype,".Rdata")
load(datdir)

datname = "case2"
dat = final_case2_imp
y = as.numeric(dat$y)
k = 10
kCV = createFolds(y, k = 10)

resdir = paste0("res/", dattype)
# #############################################################################################################
# ## GLM
# FITGLM = FITRES(dat, kCV, method="glm", upsamp=upsamp)
# png(paste0("res/roc_glm_", datname, "_", dattype, ".png"), width = 600, height = 600)
# cvROC(FITGLM$yPRED, FITGLM$yTRUE)
# dev.off()
# png(paste0(resdir,"_", datname, "_glm_rsd.png"), width = 1000, height = 600)
# binResids(FITGLM$yPRED, FITGLM$yTRUE, ylim= c(-3, 3))
# dev.off()
# png(paste0("res/vip_glm_", datname, "_", dattype, ".png"), width = 600, height = 600)
# plotVIlist(FITGLM$vIMP, p = 1)
# dev.off()
# 
# #############################################################################################################
# ## BOOSTING
# FITADA = FITRES(dat, kCV, method="ada", upsamp=upsamp)
# png(paste0("res/roc_ada_", datname, "_", dattype, ".png"), width = 600, height = 600)
# cvROC(FITADA$yPRED, FITADA$yTRUE)
# dev.off()
# png(paste0(resdir,"_", datname, "_ada_rsd.png"), width = 1000, height = 600)
# binResids(FITADA$yPRED, FITADA$yTRUE, ylim= c(-3, 3))
# dev.off()
# png(paste0("res/vip_ada_", datname, "_", dattype, ".png"), width = 600, height = 600)
# plotVIlist(FITADA$vIMP, p = 1)
# dev.off()

#############################################################################################################
## RANDOM FOREST
FITRF = FITRES(dat, kCV, method="rf", upsamp=upsamp)
png(paste0("res2/image/roc_rf_", datname, "_", dattype, ".png"), width = 600, height = 600)
cvROC(FITRF$yPRED, FITRF$yTRUE)
dev.off()
png(paste0("res2/image/acc_rf_", datname, "_", dattype, ".png"), width = 600, height = 600)
cvACC(FITRF$yPRED, FITRF$yTRUE)
dev.off()
png(paste0("res2/image/resi_rf_", datname, "_", dattype, ".png"), width = 1000, height = 600)
binResids(FITRF$yPRED, FITRF$yTRUE, ylim= c(-3, 3))
dev.off()
png(paste0("res2/image/vip_rf_", datname, "_", dattype, ".png"), width = 600, height = 600)
plotVIlist(FITRF$vIMP, p = 1)
dev.off()
#############################################################################################################
# final
final_fit = randomForest(as.factor(y)~., data = dat)
predprob = predict(final_fit, newdata = dat, type="prob")                               
cvACC(list(predprob[,2]), list(dat$y))
# acc가 너무 높은데?


# # test..
# prop = 0.8
# split = sample(1:nrow(final_case4_imp),nrow(final_case4_imp)*prop,replace=FALSE)
# train = final_case4_imp[split,]
# test = final_case4_imp[-split,]
# train_fit = randomForest(as.factor(y)~., data = train)
# test_prob = predict(train_fit, newdata = train, type="prob")                               
# cvACC(list(test_prob[,2]), list(train$y))
# 

#############################################################################################################
# save
optim_cutoff = 0.51
pred = ifelse(predprob[,2] > optim_cutoff,1,0)
true = final_case2_imp$y
mat = table(pred,true)

# label
load("data/case_data_imp_ind.Rdata")

case2 = createWorkbook(datname)
addWorksheet(case2,"predprob")
addWorksheet(case2,"predlabel")
addWorksheet(case2,"ConfusionMatrix")

# dataframe
df1 = data.frame(predprob[,2])
df1$ind = final_case2_imp_ind$ind

df2 = data.frame(data.frame(pred))
df2$ind = final_case2_imp_ind$ind


writeDataTable(case2,"predprob",df1)
writeDataTable(case2,"predlabel",df2)
writeDataTable(case2,"ConfusionMatrix",data.frame(as.matrix(mat)))

saveWorkbook(case2,file=paste0("res2/xlsx/",datname,".xlsx"),overwrite = T)

#############################################################################################################

save(FITRF,file=paste0("data/rf",datname,"_res2.Rdata"))
FITRF$yPRED


