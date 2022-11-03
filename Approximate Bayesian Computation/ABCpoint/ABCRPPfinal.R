setwd("name of directory")
library(xlsx)
library(spatstat)
library(glmnet)
###############################
#Simulation data for Pilot Run#
###############################
L=10000 #number of simulated datasets
R=100  #number of evaluation grid for K 
g=seq(0.001, 0.1, length=R)
K=seq(0.001, 0.1, length=R)
#########
#Strauss#
#########
for(k in 1:L){
  beta=runif(1,50,400)
  gamma=runif(1,0,1)
  X=rStrauss(beta,gamma,0.05,square(1))
  Kfunc=as.function(Kest(X, correction="isotropic"))
  for(r in 1:R){ 
    K[r]=Kfunc(g[r])
  }   
  X=data.frame(X)
  N=nrow(X)
  write.table(beta, "beta.txt", quote=F, 
             col.names=F, append=T)
  write.table(gamma, "gamma.txt", quote=F, 
            col.names=F, append=T)
}
###########
#DPP-Gauss#
###########
for(k in 1:L){
      tau=rgamma(1,1000,2)
      ALM=1/sqrt(pi*tau) #Gaussian
      sigma=runif(1,0.001,ALM)
      m=dppGauss(lambda=tau, alpha=sigma, d=2)
      X=simulate(m,1)
      Kfunc=as.function(Kest(X, correction="isotropic"))
      for(r in 1:R){ 
        K[r]=Kfunc(g[r])
      }   
      X=data.frame(X)
      N=nrow(X)
  write.table(tau, "tau.txt", quote=F, 
              col.names=F, append=T)
  write.table(sigma, "sigma.txt", quote=F, 
              col.names=F, append=T)
}
#######################
#DPP-Power-Exponential#
#######################
for(k in 1:L){
      tau=rgamma(1,1000,2)
      nu=10
      ALM=sqrt(gamma(2/nu+1)/(tau*gamma(2)/pi))
      alpha=runif(1,0.001,ALM)
      m=dppPowerExp(lambda=tau, alpha=alpha, nu=nu, d=2)
      X=simulate(m,1)
      Kfunc=as.function(Kest(X, correction="isotropic"))
      for(r in 1:R){ 
        K[r]=Kfunc(g[r])
      }   
      X=data.frame(X)
      N=nrow(X)
  write.table(tau, "tau.txt", quote=F, 
              col.names=F, append=T)
  write.table(alpha, "alpha.txt", quote=F, 
              col.names=F, append=T)
  write.table(N, "N.txt", quote=F, 
              col.names=F, append=T)
  write.table(X, "X.txt", quote=F, 
              col.names=F, append=T)
  write.table(K, "K.txt", quote=F, 
              col.names=F, append=T)
}


###########
#Pilot-Run#
###########
#Data format#
#  g    K 
#0.001  0
#0.002  0
#  .    .
#  .    .
#0.100 0.2
Data<-read.csv("SimF_Strauss.csv",header=T) #read the coal_ash data
Data=data.frame(Data)
N_ob=nrow(Data) #number of observed/true points
DS=ppp(Data$X,Data$Y)
R=100
g=seq(0.001, 0.1, length=R)
K=seq(0.001, 0.1, length=R)
Kfunc=as.function(Kest(DS, correction="isotropic"))
for(r in 1:R){ 
  K[r]=Kfunc(g[r])
}   
#K-function of pilot run
#Data format# (R*L, 1) vector
Sim<-read.csv("SimF_Strauss_Pilot.csv",header=T) 
Sim=data.frame(Sim)
#Paramter for pilot run 
#Data format# (L, p+1) vector
#(example: Strauss)
# beta gamma n 
# 200   0.1  90
# 198   0.5  120
#  .     .    .
#  .     .    .
SimPara<-read.csv("SimF_Strauss_Para_Pilot.csv",header=T) 
SimPara=data.frame(SimPara)
SimPara=as.matrix(SimPara)
L=nrow(SimPara) 
M=10 #number of evaluation grids, set M=1 for Strauss process 
#Observed/True Summary Statistics vector: (M+1,1) vector
Sfunc_obR=seq(0, 1, length=M+1)
for(i in 1:M){  
    Sfunc_obR[i]=sqrt(K[(R/M)*i])
}
Sfunc_obR[M+1]=N_ob

#Distance
vD=matrix(0, nrow=L, ncol=M+2)
#Pilot Summary Statistics vector: (M+1,1) vector
Sfunc_simR=seq(0, 1, length=M+1)

for(i in 1:L){
  for(j in 1:M){
    Sfunc_simR[j]=sqrt(Sim[(R*(i-1))+(R/M)*j,1])
  }
  Sfunc_simR[M+1]=SimPara[i,3]
  vD[i,1]=1
  vD[i,2:(M+1)]=(Sfunc_simR[1:M]-Sfunc_obR[1:M])^2
  vD[i,M+2]=log(Sfunc_simR[M+1])-log(Sfunc_obR[M+1])
}
betahat=solve(t(vD)%*%vD)%*%t(vD)%*%log(SimPara[,1:3])
gnet=glmnet(x=vD[,2:(M+2)],y=log(SimPara[,1:2]),family="mgaussian",alpha=1)
cv.gnet=cv.glmnet(x=vD[,2:(M+2)],y=log(SimPara[,1:2]),family="mgaussian",alpha=1)
coef(gnet,s=cv.gnet$lambda.min)
Coeff=coef(cv.gnet,s="lambda.min")
#Distance measure
vDs=seq(0, 1, length=L)
#########
#Strauss#
#########
Beta=as.matrix(Coeff$beta)
Gamma=as.matrix(Coeff$gamma)
VarBeta=var(vD[,]%*%Beta)
VarGamma=var(vD[,]%*%Gamma)
for(i in 1:L){
  vDs[i]=(Beta[1]-vD[i,]%*%Beta)^2/VarBeta+(Gamma[1]-vD[i,]%*%Gamma)^2/VarGamma
}
######
#DPPG#
######
Tau=as.matrix(Coeff$tau)
Sigma=as.matrix(Coeff$sigma)
VarTau=var(vD[,]%*%Tau)
VarSigma=var(vD[,]%*%Sigma)
for(i in 1:L){
  vDs[i]=(Tau[1]-vD[i,]%*%Tau)^2/VarTau+(Sigma[1]-vD[i,]%*%Sigma)^2/VarSigma
}
#######
#DPPPE#
#######
Tau=as.matrix(Coeff$tau)
Alpha=as.matrix(Coeff$alpha)
VarTau=var(vD[,]%*%Tau)
VarAlpha=var(vD[,]%*%Alpha)
for(i in 1:L){
  vDs[i]=(Tau[1]-vD[i,]%*%Tau)^2/VarTau+(Alpha[1]-vD[i,]%*%Alpha)^2/VarAlpha
}
vDorder=vDs[order(vDs)]
Q1=vDorder[L/100] #Set to 1% quantile 

######################
#ABC-MCMC for Strauss#
######################
Data<-read.csv("SimF_Strauss.csv",header=T) 
Data=data.frame(Data)
N_ob=nrow(Data) #the number of observed/true data points
DS=ppp(Data$X,Data$Y)
Kfunc_ob=as.function(Kest(DS, correction="isotropic"))
#Coefficient Alpha and Beta estimated through pilot run
#(M+2, p) vector
#Data format#
# beta gamma 
# 5.210 -1.780
# -291.2 657.9
#  .    .
#  .    .
theta<-read.csv("SimF_Strauss_loglasso_M10_beta.csv",header=T) 
theta=data.frame(theta)
T=1000 #the number of ABC-MCMC iteration
#Posterior Samples for Parameters
Npos=seq(0, 1, length=T+1) 
betapos=seq(0, 1, length=T+1) 
gammapos=seq(0, 1, length=T+1)
#Initial Values
betapos[1]=exp(theta[1,1])
gammapos[1]=exp(theta[1,2])
M=10 #number of functional summary statistics
g=seq(0.1/M, 0.1, length=M)
vD=seq(0, 0, length=M+2)
vD[1]=1
Kfunc_obR=seq(0, 1, length=M+1)
for(i in 1:M){
  Kfunc_obR[i]=Kfunc_ob(g[i])
}
Kfunc_obR[M+1]=log(N_ob)
for(k in 1:T){
  repeat{
    beta=betapos[k]+rnorm(1,0,20)
    gamma=runif(1,0,1)
    X=rStrauss(beta,gamma,0.05,square(1))
    Kfunc=as.function(Kest(X, correction="isotropic"))
    X=data.frame(X)
    N=nrow(X)
    Kfunc_simR=seq(0, 1, length=M+1)
    Kfunc_simR[1:M]=Kfunc(g[1:M])
    Kfunc_simR[M+1]=log(N)
    vD[1]=1
    for(j in 1:M){
     vD[j+1]=(Kfunc_simR[j]^(1/2)-Kfunc_obR[j]^(1/2))^2
    }
     vD[M+2]=(Kfunc_simR[M+1]-Kfunc_obR[M+1])
    EZ=(theta[1,1]-t(theta[,1])%*%vD)^2/VarBeta+(theta[1,2]-t(theta[,2])%*%vD)^2/VarGamma
    if(EZ<=Q1){
      betapos[k+1]=beta
      gammapos[k+1]=gamma
      Npos[k+1]=N
      break
    }
  }
  write.table(betapos[k+1], "betapos.txt", quote=F, 
              col.names=F, append=T)
  write.table(gammapos[k+1], "gammapos.txt", quote=F, 
              col.names=F, append=T)
  write.table(Npos[k+1], "Npos.txt", quote=F, 
              col.names=F, append=T)
  write.table(X, "X.txt", quote=F, 
            col.names=F, append=T)
}

###############
#ABC-MCMC-DPPG#
###############
#location of observed/true points
Data<-read.csv("SimF_DPP.csv",header=T) 
Data=data.frame(Data)
N_ob=nrow(Data) #the number of observed/true data points
DS=ppp(Data$X,Data$Y)
Kfunc_ob=as.function(Kest(DS, correction="isotropic"))
#Coefficient Alpha and Beta estimated through pilot run
#(M+2, p) vector
#Data format#
#  tau   sigma 
# 5.210 -1.780
# -291.2 657.9
#  .    .
#  .    .
theta<-read.csv("SimF_DPP_loglasso_M10_beta.csv",header=T) #read the coal_ash data
theta=data.frame(theta)
T=1000
Npos=seq(0, 1, length=T+1)
taupos=seq(0, 1, length=T+1)
sigmapos=seq(0, 1, length=T+1)
tau=theta[1,1]
sigma=theta[1,2]
M=10 #number of evaluation grids for functional summary statistics
g=seq(0.1/M, 0.1, length=M)
Kfunc_obR=seq(0, 1, length=M+1)
for(i in 1:M){
  Kfunc_obR[i]=Kfunc_ob(g[i])
}
Kfunc_obR[M+1]=log(N_ob)
vD=seq(0, 0, length=M+2)
vD[1]=1
for(k in 1:T){
  repeat{
    tau=taupos[k]+rnorm(1,0,5)
    ALM=1/sqrt(pi*tau) #Gaussian
    sigma=runif(1,0.001,ALM)
    m=dppGauss(lambda=tau, alpha=sigma, d=2)
    X=simulate(m,1)
    X=data.frame(X)
    N=nrow(X)
    Kfunc_simR=seq(0, 1, length=M+1)
    Kfunc_simR[1:M]=Kfunc(g[1:M])
    Kfunc_simR[M+1]=log(N)
    vD[1]=1
    for(j in 1:M){
      vD[j+1]=(Kfunc_simR[j]^(1/2)-Kfunc_obR[j]^(1/2))^2
    }
    vD[M+2]=(Kfunc_simR[M+1]-Kfunc_obR[M+1])
    EZ=(theta[1,1]-t(theta[,1])%*%vD)^2/VarTau+(theta[1,2]-t(theta[,2])%*%vD)^2/VarSigma
    if(EZ<=Q1){
      taupos[k+1]=tau
      sigmapos[k+1]=sigma
      Npos[k+1]=N
      break
    }
  }
  write.table(taupos[k+1], "taupos.txt", quote=F, 
              col.names=F, append=T)
  write.table(sigmapos[k+1], "sigmapos.txt", quote=F, 
              col.names=F, append=T)
  write.table(Npos[k+1], "Npos.txt", quote=F, 
              col.names=F, append=T)
  write.table(X, "X.txt", quote=F, 
              col.names=F, append=T)
}

################
#ABC-MCMC-DPPPE#
################
#location of observed/true points
Data<-read.csv("SimF_DPPPE.csv",header=T) 
Data=data.frame(Data)
N_ob=nrow(Data) #the number of observed/true data points
DS=ppp(Data$X,Data$Y)
Kfunc_ob=as.function(Kest(DS, correction="isotropic")) #K-function for observed/true points
#Coefficient Alpha and Beta estimated through pilot run
#(M+2, p) vector
theta<-read.csv("SimF_DPPPE_loglasso_M10_beta.csv",header=T) 
theta=data.frame(theta)
M=10 #number of evaluation grids for functional summary statistics
T=1000
Npos=seq(0, 1, length=T+1)
taupos=seq(0, 1, length=T+1)
alphapos=seq(0, 1, length=T+1)
tau=theta[1,1]
alpha=theta[1,2]
nu=10 #parameter for DPP-PE
g=seq(0.1/M, 0.1, length=M)
Kfunc_obR=seq(0, 1, length=M+1)
for(i in 1:M){
  Kfunc_obR[i]=Kfunc_ob(g[i])
}
Kfunc_obR[M+1]=log(N_ob)
vD=seq(0, 0, length=M+2)
for(k in 1:T){
  repeat{
    tau=taupos[k]+rnorm(1,0,5)
    ALM=sqrt(gamma(2/nu+1)/(tau*gamma(2)/pi))
    alpha=runif(1,0.001,ALM)
    m=dppPowerExp(lambda=tau, alpha=alpha, nu=nu, d=2)
    X=simulate(m,1)
    X=data.frame(X)
    N=nrow(X)
    Kfunc_simR=seq(0, 1, length=M+1)
    Kfunc_simR[1:M]=Kfunc(g[1:M])
    Kfunc_simR[M+1]=log(N)
    vD[1]=1
    for(j in 1:M){
      vD[j+1]=(Kfunc_simR[j]^(1/2)-Kfunc_obR[j]^(1/2))^2
    }
    vD[M+2]=(Kfunc_simR[M+1]-Kfunc_obR[M+1])
    EZ=(theta[1,1]-t(theta[,1])%*%vD)^2/VarTau+(theta[1,2]-t(theta[,2])%*%vD)^2/VarAlpha
    if(EZ<=Q1){
      taupos[k+1]=tau
      alphapos[k+1]=alpha
      Npos[k+1]=N
      break
    }
  }
  write.table(taupos[k+1], "taupos.txt", quote=F, 
              col.names=F, append=T)
  write.table(alphapos[k+1], "alphapos.txt", quote=F, 
              col.names=F, append=T)
  write.table(Npos[k+1], "Npos.txt", quote=F, 
              col.names=F, append=T)
  write.table(X, "X.txt", quote=F, 
              col.names=F, append=T)
}
