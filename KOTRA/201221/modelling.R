rm(list=ls())
setwd("C:/Users/user/Seungjun/KOTRA/201221")

source("helper.r")
library(tidyverse)
library(openxlsx)
#############################################################################################################
## create k-cv folds
dattype = "imp"
upsamp=F
datdir = paste0("data/", "case_data_", dattype,".Rdata")
load(datdir)

datname = "case4"
dat = final_case4_imp
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
png(paste0("res/image/roc_rf_", datname, "_", dattype, ".png"), width = 600, height = 600)
cvROC(FITRF$yPRED, FITRF$yTRUE)
dev.off()
png(paste0("res/image/acc_rf_", datname, "_", dattype, ".png"), width = 600, height = 600)
cvACC(FITRF$yPRED, FITRF$yTRUE)
dev.off()
png(paste0("res/image/resi_rf_", datname, "_", dattype, ".png"), width = 1000, height = 600)
binResids(FITRF$yPRED, FITRF$yTRUE, ylim= c(-3, 3))
dev.off()
png(paste0("res/image/vip_rf_", datname, "_", dattype, ".png"), width = 600, height = 600)
plotVIlist(FITRF$vIMP, p = 1)
dev.off()
#############################################################################################################
# save

mat = ConfusionMat(FITRF$yPRED, FITRF$yTRUE,0.51)
final

write.xlsx(mat,sheetName=predprob,file="res/xlsx/result.xlsx")
write.xlsx(mat,sheetName=predlabel,file="res/xlsx/result.xlsx")
write.xlsx(mat,sheetName=ConfusionMatrix,file="res/xlsx/result.xlsx")
#############################################################################################################

save(FITRF,file="data/rfcase4_res.Rdata")
FITRF$yPRED


