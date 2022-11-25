library(sp)
library(gstat)
library(fields)
library(classInt)
library(spatstat)
library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(spatstat)
library(fields)
library(mvtnorm)
setwd("C:/Users/user/Seungjun/GCDS/Simulation")
#-----------------------------------------------------------------------------#
n = 30
x = rep(0:(n - 1) / (n - 1), times = n)
y = rep(0:(n - 1) / (n - 1), each = n)
X = cbind(x, y)
dist_X = rdist(X)

t = 30
t_interval = 1:t
t_interval =t_interval/length(t_interval)
dist_t = rdist(t_interval)

# spatial correlation
phi_spat = 1 #  # positive real number
delta_spat = 1 # (0,2]
rho_spat = exp(-(dist_X/phi_spat)^(2+delta_spat))
dim(rho_spat)
# temporal correlation
phi_temp = 0.5 #  [-1,1]
gamma_temp = phi_temp^(dist_t)
dim(gamma_temp)


# spatio temporal correlation
st_corr = base::kronecker(gamma_temp,rho_spat)
dim(st_corr)

#-----------------------------------------------------------------------------#

set.seed(10)
eta = mgcv::rmvn(1, rep(0,n*n*t), 2*st_corr)
eta_mat = matrix(eta,nrow=t,ncol=n*n,byrow=TRUE)


par(mfrow=c(5,6),mar=c(1,0,1.5,0))
for(i in 1:t){
  image(matrix(eta_mat[i,],ncol=n,byrow=TRUE),main=paste0("corr;t",i),axes = FALSE)
  
}

# eta = matrix(sample,ncol=900,byrow=TRUE)
length(eta)

X = mgcv::rmvn(n*n*t,rep(0,5),0.005*diag(5))
beta = c(1,1,1,1,1)
Xbeta = X%*%beta
dim(Xbeta)
y = Xbeta + eta

st_samples = matrix(y,nrow=t,ncol=n*n,byrow=TRUE)
dim(st_samples)


par(mfrow=c(5,6),mar=c(1,0,1.5,0))
for(i in 1:t){
  image(matrix(st_samples[i,],ncol=n,byrow=TRUE),axes=FALSE,main=paste0("y;t",i))
}

#-----------------------------------------------------------------------------#
# lag
len = dim(st_samples)[1]
idx_lst = combn(1:len,2)
lag = idx_lst[2,] - idx_lst[1,]
lag


#-----------------------------------------------------------------------------#
# making pairs
# pass

write.csv(st_samples,"Simulation_data.csv",row.names=FALSE)
?write.csv





