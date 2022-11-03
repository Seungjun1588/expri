setwd("C:/Users/user/Desktop/승준/연구/Sx_DNN")
library(Rcpp)
library(RcppArmadillo)
sourceCpp("Ising_abc.cpp")

nx =10; ny=10
# theta =0.3을 기준으로 생성
# Obs -----------------------------------------------------------------------#
iter = 30

set.seed(444)
obs = Gen_data(size=c(nx,ny),n=iter)
write.csv(obs$X,"obs.csv",row.names=FALSE)
write.csv(obs$theta,"obs_theta.csv",row.names=FALSE)
rm(obs)

# train ---------------------------------------------------------------------#
iter = 1e+06

set.seed(111)
train_data = Gen_data(size=c(nx,ny),n=iter)
write.csv(train_data$X,"train/train_x.csv",row.names=FALSE)
write.csv(train_data$theta,"train/train_theta.csv",row.names=FALSE)
write.csv(train$Summary_stat,"train/train_SS.csv",row.names=FALSE)
rm(train_data)

# val -----------------------------------------------------------------------#
iter = 1e+05

set.seed(222)
val_data = Gen_data(size=c(nx,ny),n=iter)
write.csv(val_data$X,"val/val_x.csv",row.names=FALSE)
write.csv(val_data$theta,"val/val_theta.csv",row.names=FALSE)
write.csv(val_data$Summary_stat,"val/val_SS.csv",row.names=FALSE)
rm(val_data)

# test ----------------------------------------------------------------------#
iter = 1e+05

set.seed(333)
test_data = Gen_data(size=c(nx,ny),n=iter)
write.csv(test_data$X,"test/test_x.csv",row.names=FALSE)
write.csv(test_data$theta,"test/test_theta.csv",row.names=FALSE)
write.csv(test_data$Summary_stat,"test/test_SS.csv",row.names=FALSE)
rm(test_data)
