#########--------------------------------------------------------------#######
#########--------------------------------------------------------------#######
#########--------------------------------------------------------------#######
#########--------------------------------------------------------------#######

setwd("C:/Users/user/desktop/jun/SpatioTemp/lecture4")

source("http://www.stat.psu.edu/~mharan/batchmeans.R")
library(nimble);library(mvtnorm);library(fields)
library(classInt)
library(fields)
library(maps)
library(sp)
library(gstat)
library(geoR)
library(mvtnorm)
library(MCMCpack)
library(coda)

ploteqc <- function(spobj, z, breaks, ...){
  pal <- tim.colors(length(breaks)-1)
  fb <- classIntervals(z, n = length(pal), 
                       style = "fixed", fixedBreaks = breaks)
  col <- findColours(fb, pal)
  plot(spobj, col = col, ...)
  image.plot(legend.only = TRUE, zlim = range(breaks), col = pal)
}

## Using Ming-Hui Chen's paper in Journal of Computational and Graphical Stats.
hpd <- function(samp,p=0.05){
  ## to find an approximate (1-p)*100% HPD interval from a
  ## given posterior sample vector samp
  
  r <- length(samp)
  samp <- sort(samp)
  rang <- matrix(0,nrow=trunc(p*r),ncol=3)
  dimnames(rang) <- list(NULL,c("low","high","range"))
  for (i in 1:trunc(p*r)) {
    rang[i,1] <- samp[i]
    rang[i,2] <- samp[i+(1-p)*r]
    rang[i,3] <- rang[i,2]-rang[i,1]
  }
  hpd <- rang[order(rang[,3])[1],1:2]
  return(hpd)
}

# Exponential Covariance Function
expCov<-function(distMat,phi){
  exp(-distMat/phi)
}

sqeCov<-function(distMat,phi){
  exp(-0.5*(distMat/phi)^2)
}

matCov<-function(distMat,phi){
  (1+(sqrt(5)*(distMat/phi))+((5*distMat^2)/(3*(phi^2))))*exp(-(sqrt(5)*(distMat/phi)))
}


# Matern Cov Function + Acceptance Rate function
Matern <- function(d, param = c(scale = 1, range = 1, smoothness = 2)) {
  scale <- param[1]
  range <- param[2]
  smoothness <- param[3]
  if (any(d < 0))
    stop("distance argument must be nonnegative")
  d <- d / range
  d[d == 0] <- 1e-10
  rootcon<-sqrt(2*smoothness)
  con <- (2^(smoothness - 1)) * gamma(smoothness)
  con <- 1 / con
  return(scale * con * ((rootcon*d)^smoothness) * besselK(rootcon*d, smoothness))
}
accRateFunc<-function(x){
  accRate<-(length(unique(x))-1)/(length(x)-1)
  return(accRate)
}


# Summary 
summaryFunction<-function(mcmcDat,bmseThresh=0.01,time){
  
  # Parameters
  summaryMat<-rbind(apply(mcmcDat,2,mean),
                    apply(mcmcDat,2,hpd),
                    apply(mcmcDat,2,accRateFunc),
                    bmmat(mcmcDat)[,2],
                    abs(apply(mcmcDat,2,mean))*bmseThresh,
                    apply(mcmcDat,2,ess),
                    apply(mcmcDat,2,ess)/time)
  
  rownames(summaryMat)<-c("Mean","95%CI-Low","95%CI-High",
                          "Accept","BMSE",paste(bmseThresh,"x mean"),
                          "ESS","ESS/sec")
  return(summaryMat)
}

# NIMBLE FUNCTIONS

expcov <- nimbleFunction(     
  run = function(dists = double(2), phi = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    for(i in 1:n){
      for(j in 1:n){
        result[i, j] <- exp(-dists[i,j]/phi)
      }
    }
    
    return(result)
  })


matcov <- nimbleFunction(     
  run = function(dists = double(2), phi = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    for(i in 1:n){
      for(j in 1:n){
        result[i, j] <- (1+(sqrt(5)*(dists[i,j]/phi))+((5*dists[i,j]^2)/(3*(phi^2))))*exp(-(sqrt(5)*(dists[i,j]/phi)))
      }
    }
    
    return(result)
  })


#########--------------------------------------------------------------#######
#########--------------------------------------------------------------#######
#########--------------------------------------------------------------#######
#########--------------------------------------------------------------#######
# (a)

######### 2. Simulate Dataset ###############

#Parameters

Simulate_binary_data = function(n=400,cv=0.2){
  set.seed(55555)
  ncv=n*cv
  beta=c(1,1) ; phi=0.2 ; sigma2=1
  
  # Generate Locations 
  ## Split into Model + Cross-Validation
  gridLocation<-cbind(runif(n,min = 0,max = 1),runif(n,min = 0,max = 1)) # 400*2
  CVgridLocation<-cbind(runif(ncv,min = 0,max = 1),runif(ncv,min = 0,max = 1)) # 100*2
  comboLocation<-rbind(gridLocation,CVgridLocation)  # 500*2
  distMatFull<-as.matrix(rdist(comboLocation)) # 500*500
  # Create Indices
  modInd<-1:n # 1~400
  CVInd<-(n+1):nrow(distMatFull) # 401~500
  distMatMod<-distMatFull[modInd,modInd] # 400*400
  # Covariates 
  XMat<-cbind(runif(n,-1,1),runif(n,-1,1)) # 400*2
  XMatCV<-cbind(runif(ncv,-1,1),runif(ncv,-1,1)) # 100*2
  XB<-XMat%*%beta # 400*1
  cvXB<-XMatCV%*%beta #100*1
  XBFull<-rbind(XB,cvXB) # 500*1
  
  # Covariance Matrix
  CovMat<-sigma2*matCov(distMatFull,phi)
  
  # Latent Gaussian Random Field
  gpWFull <- as.numeric(rmvnorm(n=1,mean=rep(0,nrow(CovMat)),sigma = CovMat,method = "chol"))
  pWFullLinear<-gpWFull+XBFull
  pWFullBin<-exp(gpWFull+XBFull)/(1+exp(gpWFull+XBFull)) # 500*1
  
  # Observations
  obsFullLinear<-pWFullLinear # 500*1
  obsFullBin<-sapply(pWFullBin,rbinom,n=1,size=1)
  
  
  ##################
  # Model Sample
  # Latent Process
  gpWMod<-gpWFull[modInd]
  # Expected Value
  pWModLinear<-pWFullLinear[modInd] #400*1
  pWModBin<-pWFullBin[modInd]
  
  # Observations
  obsModLinear<-obsFullLinear[modInd]
  obsModBin<-obsFullBin[modInd]
  
  # CV Sample
  # Latent Process
  gpWCV<-gpWFull[CVInd]
  # Expected Value
  pWCVLinear<-pWFullLinear[CVInd]
  pWCVBin<-pWFullBin[CVInd]
  # Observations
  obsCVLinear<-obsFullLinear[CVInd] # 100*1
  obsCVBin<-obsFullBin[CVInd]
  
  return(list("obsModBin"=obsModBin,"obsCVBin"=obsCVBin,"XMat"=XMat,"XmatCV"=XMatCV,
              "distMatMod"=distMatMod,"gridLocation"=gridLocation,"CVgridLocation"=CVgridLocation))
}

#########--------------------------------------------------------------#######
#########--------------------------------------------------------------#######
#########--------------------------------------------------------------#######
#########--------------------------------------------------------------#######
# (b)


### for n= 200

model_string <- nimbleCode({
  
  # Data Model
  for(i in 1:n){
    lambda[i] <- exp(W[i]+XB[i])
    Z[i] ~ dpois(lambda[i])
  }
  
  # Constant and Cov Matrix
  XB[1:n]<-beta1*X[,1] + beta2*X[,2]
  covMat[1:n,1:n]<- expcov(dists[1:n,1:n],phi)
  fullCovMat[1:n,1:n]<- sigma2*covMat[1:n,1:n]
  
  # Process Model
  W[1:n] ~ dmnorm(mean = mn[1:n], cov = fullCovMat[1:n,1:n])
  
  # Parameter Model
  sigma2   ~  dinvgamma(0.2, 0.2)
  phi   ~  dunif(0,1)
  beta1 ~  dnorm(0, sd=sqrt(100))
  beta2 ~  dnorm(0, sd=sqrt(100))
})

niter=50000
n =200
sim.data = Simulate_binary_data(n=n,cv=0.2)
obsModBinom = sim.data$obsModBin
XMat = sim.data$XMat
distMatMod = sim.data$distMatMod

consts   <- list(n=n,X=XMat,dists=distMatMod,mn=rep(0,n))
data     <- list(Z=obsModBinom)
inits    <- list(beta1=rnorm(1),beta2=rnorm(1),phi=0.5,sigma2=2, 
                 W=rnorm(n))

# Run MCMC


nimble_model = nimbleModel(model_string, data = data, inits = inits,constants=consts)
MCMC_conf = configureMCMC(nimble_model,monitors = c("beta1", "beta2","phi","sigma2"),
                          control = list(adaptFactorExponent = 1))
RMCMC = buildMCMC(MCMC_conf,samplesAsCodaMCMC=TRUE,WAIC=FALSE,summary=FALSE, thin=20,
                  niter = niter, nburnin = 0, nchains = 1)
Cmodel = compileNimble(nimble_model)
Cmcmc = compileNimble(RMCMC,project = nimble_model)

pt<-proc.time()
Cmcmc$run(niter = niter)
samples200 = as.matrix(Cmcmc$mvSamples)

# samples200  <- nimbleMCMC(model_string, data = data, inits = inits,
#                           constants=consts,
#                           monitors = c("beta1", "beta2","phi","sigma2"),
#                           samplesAsCodaMCMC=TRUE,WAIC=FALSE,summary=FALSE, thin=20,
#                           niter = niter, nburnin = 0, nchains = 1)
ptFinal200<-proc.time()-pt

# computing time
ptFinal200

# trace plot
pdf(file = "BinomResults200.pdf",width=11,height=8.5)
par(mfrow=c(4,2),mar=c(2,2,2,2))
sampInd<-floor(seq(1,nrow(samples200),length.out = 1000))
beta=c(1,1) ; phi=0.2 ; sigma2=1

for(i in 1:4){
  plot.ts(samples200[sampInd,i]); abline(h=c(beta,phi,sigma2)[i],col="red",lwd=2)
  plot(density(samples200[sampInd,i])); abline(v=c(beta,phi,sigma2)[i],col="red",lwd=2)
}

dev.off()

summaryMat200<-list()
summaryMat200[[1]]<-round(summaryFunction(samples200[,c("beta1", "beta2","phi","sigma2")],
                                          time=ptFinal200[3]),3)
# posterior mean, hpd, acc rate
summaryMat200[[1]]
save(samples200,file="BinomMCMCsamples200.RData")
save(summaryMat200,samples200,ptFinal200,file="BinomMCMCResults200.RData")

# prediction acc
burnin <- 100
s2.final <- samples200[,c("sigma2")][-(1:burnin)]
beta.final <- cbind(samples200[,c("beta1")][-(1:burnin)],samples200[,c("beta2")][-(1:burnin)])
rho.final <- samples200[,c("phi")][-(1:burnin)]


obs.grid = sim.data$gridLocation
pred.grid = sim.data$CVgridLocation
X = sim.data$XMat
Xpred = sim.data$XmatCV

full.grid = rbind(obs.grid,pred.grid)
full.X = rbind(X,Xpred)
# fixed part
full.XB = full.X%*%t(beta.final)

eta.mat <- matrix(NA, nrow = nrow(pred.grid), ncol = niter-burnin) # 40*

# random part
distMatFull = rdist.earth(coordinates(full.grid))

for(j in 1:ncol(eta.mat)){
  if(j%%200 ==0){print(j)}
  
  CovMat.full <-s2.final[j]*matCov(distMatFull,rho.final[j])
  gpWFull <- as.numeric(rmvnorm(n=1,mean=rep(0,nrow(CovMat.full)),sigma = CovMat.full,method = "chol"))
  eta.mat[,j] <- (gpWFull+full.XB[,j])[(n+1):length(gpWFull)]
  
}
save(eta.mat,file="Binom_eta200.RData")
dim(eta.mat)

acc_lst = NULL
actual.y = sim.data$obsCVBin

for(i in 1:ncol(eta.mat)){
  bin.p = exp(eta.mat[,i])/(1+exp(eta.mat[,i]))
  pred.y = sapply(bin.p,rbinom,n=1,size=1)
  ct = table(pred.y,actual.y)
  acc_lst[i] = sum(diag(ct))/sum(ct)

}
load("Binom_eta200.RData")
save(acc_lst,file="acc_lst200.RData")

length(acc_lst)
mean(acc_lst)


