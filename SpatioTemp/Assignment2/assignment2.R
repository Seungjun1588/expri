#########--------------------------------------------------------------#######
#########--------------------------------------------------------------#######
#########--------------------------------------------------------------#######
#########--------------------------------------------------------------#######

setwd("C:/Users/user/desktop/jun/SpatioTemp/lecture4")

source("http://www.stat.psu.edu/~mharan/batchmeans.R")
library(nimble);library(mvtnorm);library(fields)

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

Simulate_binary_data = function(n=200,cv=0.2){
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

niter=100000
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

pt<-proc.time()
samples200  <- nimbleMCMC(model_string, data = data, inits = inits,
                       constants=consts,
                       monitors = c("beta1", "beta2","phi","sigma2"),
                       samplesAsCodaMCMC=TRUE,WAIC=FALSE,summary=FALSE,
                       niter = niter, nburnin = 0, nchains = 1)
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
summaryMat200[[2]]<-round(summaryFunction(samples200[,1:n],
                                       time=ptFinal200[3]),3)
# psterior mean, hpd, acc rate
summaryMat200[[1]]
apply(summaryMat200[[2]],1,mean)
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
y = sim.data$obsCVBin
  
d <- rdist.earth(coordinates(obs.grid))
dcross <- rdist.earth(coordinates(obs.grid), coordinates(pred.grid)) # 200*40
dpred <- rdist.earth(coordinates(pred.grid)) # 40*40

eta.pred <- matrix(NA, nrow = nrow(pred.grid), ncol = niter-burnin) # 40*

for(j in 1:ncol(eta.pred)){
  if(j%%200 == 0){print(j)}
  
  # Construct the covariance matrices
  Gamma <- exp(-d/rho.final[j]) 
  Ginv <- solve(Gamma)
  g <- exp(-dcross/rho.final[j])
  Gpred <- exp(-dpred/rho.final[j])
  
  m <- Xpred %*% beta.final[j,] + t(g) %*% Ginv %*% 
    (y - X %*% beta.final[j,])
  V <- s2.final[j] * (Gpred - t(g)%*%Ginv%*%g)
  eta.pred[,j] <- rmvnorm(1, m, V, method = "svd")
}

## Find pointwise posterior means and sds

eta.pred.m <- apply(eta.pred, 1, mean)
eta.pred.sd <- apply(eta.pred, 1, sd)

bin.p = exp(eta.pred.m)/(1+exp(eta.pred.m))

pred.y = sapply(bin.p,rbinom,n=1,size=1)
actual.y = sim.data$obsCVBin

# contingency table
ct = table(pred.y,actual.y)

# accuracy 
sum(diag(ct))/sum(ct)


#########--------------------------------------------------------------#######
#########--------------------------------------------------------------#######
#########--------------------------------------------------------------#######
#########--------------------------------------------------------------#######
# (c)
#########--------------------------------------------------------------#######

# for n = 400
 
n =400
sim.data = Simulate_binary_data(n=n,cv=0.2)
obsModBinom = sim.data$obsModBin
XMat = sim.data$XMat
distMatMod = sim.data$distMatMod

consts   <- list(n=n,X=XMat,dists=distMatMod,mn=rep(0,n))
data     <- list(Z=obsModBinom)
inits    <- list(beta1=rnorm(1),beta2=rnorm(1),phi=0.5,sigma2=2, 
                 W=rnorm(n))

# Run MCMC

pt<-proc.time()
samples400  <- nimbleMCMC(model_string, data = data, inits = inits,
                          constants=consts,
                          monitors = c("beta1", "beta2","phi","sigma2"),
                          samplesAsCodaMCMC=TRUE,WAIC=FALSE,summary=FALSE,
                          niter = niter, nburnin = 0, nchains = 1)
ptFinal400<-proc.time()-pt

# computing time
ptFinal400

# trace plot
pdf(file = "BinomResults400.pdf",width=11,height=8.5)
par(mfrow=c(4,2),mar=c(2,2,2,2))
sampInd<-floor(seq(1,nrow(samples400),length.out = 1000))
beta=c(1,1) ; phi=0.2 ; sigma2=1

for(i in 1:4){
  plot.ts(samples400[sampInd,i]); abline(h=c(beta,phi,sigma2)[i],col="red",lwd=2)
  plot(density(samples400[sampInd,i])); abline(v=c(beta,phi,sigma2)[i],col="red",lwd=2)
}

dev.off()

summaryMat400<-list()
summaryMat400[[1]]<-round(summaryFunction(samples400[,c("beta1", "beta2","phi","sigma2")],
                                          time=ptFinal400[3]),3)
summaryMat400[[2]]<-round(summaryFunction(samples400[,1:n],
                                          time=ptFinal400[3]),3)
# psterior mean, hpd, acc rate
summaryMat400[[1]]
apply(summaryMat400[[2]],1,mean)
save(samples400,file="BinomMCMCsamples400.RData")
save(summaryMat400,samples400,ptFinal400,file="BinomMCMCResults400.RData")

# prediction acc
burnin <- 100
s2.final <- samples400[,c("sigma2")][-(1:burnin)]
beta.final <- cbind(samples400[,c("beta1")][-(1:burnin)],samples400[,c("beta2")][-(1:burnin)])
rho.final <- samples400[,c("phi")][-(1:burnin)]


obs.grid = sim.data$gridLocation
pred.grid = sim.data$CVgridLocation
X = sim.data$XMat
Xpred = sim.data$XmatCV
y = sim.data$obsCVBin

d <- rdist.earth(coordinates(obs.grid))
dcross <- rdist.earth(coordinates(obs.grid), coordinates(pred.grid)) # 400*80
dpred <- rdist.earth(coordinates(pred.grid)) # 80*80

eta.pred <- matrix(NA, nrow = nrow(pred.grid), ncol = niter-burnin) # 80*400

for(j in 1:ncol(eta.pred)){
  print(j)
  
  # Construct the covariance matrices
  Gamma <- exp(-d/rho.final[j]) 
  Ginv <- solve(Gamma)
  g <- exp(-dcross/rho.final[j])
  Gpred <- exp(-dpred/rho.final[j])
  
  m <- Xpred %*% beta.final[j,] + t(g) %*% Ginv %*% 
    (y - X %*% beta.final[j,])
  V <- s2.final[j] * (Gpred - t(g)%*%Ginv%*%g)
  eta.pred[,j] <- rmvnorm(1, m, V, method = "svd")
}

## Find pointwise posterior means and sds

eta.pred.m <- apply(eta.pred, 1, mean)
eta.pred.sd <- apply(eta.pred, 1, sd)

bin.p = exp(eta.pred.m)/(1+exp(eta.pred.m))

pred.y = sapply(bin.p,rbinom,n=1,size=1)
actual.y = sim.data$obsCVBin

# contingency table
ct = table(pred.y,actual.y)

# accuracy 
sum(diag(ct))/sum(ct)




#########--------------------------------------------------------------#######

# for n = 800

n =800
sim.data = Simulate_binary_data(n=n,cv=0.2)
obsModBinom = sim.data$obsModBin
XMat = sim.data$XMat
distMatMod = sim.data$distMatMod

consts   <- list(n=n,X=XMat,dists=distMatMod,mn=rep(0,n))
data     <- list(Z=obsModBinom)
inits    <- list(beta1=rnorm(1),beta2=rnorm(1),phi=0.5,sigma2=2, 
                 W=rnorm(n))

# Run MCMC

pt<-proc.time()
samples800  <- nimbleMCMC(model_string, data = data, inits = inits,
                          constants=consts,
                          monitors = c("beta1", "beta2","phi","sigma2"),
                          samplesAsCodaMCMC=TRUE,WAIC=FALSE,summary=FALSE,
                          niter = niter, nburnin = 0, nchains = 1)
ptFinal800<-proc.time()-pt

# computing time
ptFinal800

# trace plot
pdf(file = "BinomResults800.pdf",width=11,height=8.5)
par(mfrow=c(4,2),mar=c(2,2,2,2))
sampInd<-floor(seq(1,nrow(samples800),length.out = 1000))
beta=c(1,1) ; phi=0.2 ; sigma2=1

for(i in 1:4){
  plot.ts(samples800[sampInd,i]); abline(h=c(beta,phi,sigma2)[i],col="red",lwd=2)
  plot(density(samples800[sampInd,i])); abline(v=c(beta,phi,sigma2)[i],col="red",lwd=2)
}

dev.off()

summaryMat800<-list()
summaryMat800[[1]]<-round(summaryFunction(samples800[,c("beta1", "beta2","phi","sigma2")],
                                          time=ptFinal800[3]),3)
summaryMat800[[2]]<-round(summaryFunction(samples800[,1:n],
                                          time=ptFinal800[3]),3)
# psterior mean, hpd, acc rate
summaryMat800[[1]]
apply(summaryMat800[[2]],1,mean)
save(samples800,file="BinomMCMCsamples800.RData")
save(summaryMat800,samples800,ptFinal800,file="BinomMCMCResults800.RData")

# prediction acc
burnin <- 100
s2.final <- samples800[,c("sigma2")][-(1:burnin)]
beta.final <- cbind(samples800[,c("beta1")][-(1:burnin)],samples800[,c("beta2")][-(1:burnin)])
rho.final <- samples800[,c("phi")][-(1:burnin)]


obs.grid = sim.data$gridLocation
pred.grid = sim.data$CVgridLocation
X = sim.data$XMat
Xpred = sim.data$XmatCV
y = sim.data$obsCVBin

d <- rdist.earth(coordinates(obs.grid))
dcross <- rdist.earth(coordinates(obs.grid), coordinates(pred.grid)) # 800*80
dpred <- rdist.earth(coordinates(pred.grid)) # 80*80

eta.pred <- matrix(NA, nrow = nrow(pred.grid), ncol = niter-burnin) # 80*800

for(j in 1:ncol(eta.pred)){
  print(j)
  
  # Construct the covariance matrices
  Gamma <- exp(-d/rho.final[j]) 
  Ginv <- solve(Gamma)
  g <- exp(-dcross/rho.final[j])
  Gpred <- exp(-dpred/rho.final[j])
  
  m <- Xpred %*% beta.final[j,] + t(g) %*% Ginv %*% 
    (y - X %*% beta.final[j,])
  V <- s2.final[j] * (Gpred - t(g)%*%Ginv%*%g)
  eta.pred[,j] <- rmvnorm(1, m, V, method = "svd")
}

## Find pointwise posterior means and sds

eta.pred.m <- apply(eta.pred, 1, mean)
eta.pred.sd <- apply(eta.pred, 1, sd)

bin.p = exp(eta.pred.m)/(1+exp(eta.pred.m))

pred.y = sapply(bin.p,rbinom,n=1,size=1)
actual.y = sim.data$obsCVBin

# contingency table
ct = table(pred.y,actual.y)

# accuracy 
sum(diag(ct))/sum(ct)

#########--------------------------------------------------------------#######

# for n = 1600

n =1600
sim.data = Simulate_binary_data(n=n,cv=0.2)
obsModBinom = sim.data$obsModBin
XMat = sim.data$XMat
distMatMod = sim.data$distMatMod

consts   <- list(n=n,X=XMat,dists=distMatMod,mn=rep(0,n))
data     <- list(Z=obsModBinom)
inits    <- list(beta1=rnorm(1),beta2=rnorm(1),phi=0.5,sigma2=2, 
                 W=rnorm(n))

# Run MCMC

pt<-proc.time()
samples1600  <- nimbleMCMC(model_string, data = data, inits = inits,
                          constants=consts,
                          monitors = c("beta1", "beta2","phi","sigma2"),
                          samplesAsCodaMCMC=TRUE,WAIC=FALSE,summary=FALSE,
                          niter = niter, nburnin = 0, nchains = 1)
ptFinal1600<-proc.time()-pt

# computing time
ptFinal1600

# trace plot
pdf(file = "BinomResults1600.pdf",width=11,height=8.5)
par(mfrow=c(4,2),mar=c(2,2,2,2))
sampInd<-floor(seq(1,nrow(samples1600),length.out = 1000))
beta=c(1,1) ; phi=0.2 ; sigma2=1

for(i in 1:4){
  plot.ts(samples1600[sampInd,i]); abline(h=c(beta,phi,sigma2)[i],col="red",lwd=2)
  plot(density(samples1600[sampInd,i])); abline(v=c(beta,phi,sigma2)[i],col="red",lwd=2)
}

dev.off()

summaryMat1600<-list()
summaryMat1600[[1]]<-round(summaryFunction(samples1600[,c("beta1", "beta2","phi","sigma2")],
                                          time=ptFinal1600[3]),3)
summaryMat1600[[2]]<-round(summaryFunction(samples1600[,1:n],
                                          time=ptFinal1600[3]),3)
# psterior mean, hpd, acc rate
summaryMat1600[[1]]
apply(summaryMat1600[[2]],1,mean)
save(samples1600,file="BinomMCMCsamples1600.RData")
save(summaryMat1600,samples1600,ptFinal1600,file="BinomMCMCResults1600.RData")

# prediction acc
burnin <- 100
s2.final <- samples1600[,c("sigma2")][-(1:burnin)]
beta.final <- cbind(samples1600[,c("beta1")][-(1:burnin)],samples1600[,c("beta2")][-(1:burnin)])
rho.final <- samples1600[,c("phi")][-(1:burnin)]


obs.grid = sim.data$gridLocation
pred.grid = sim.data$CVgridLocation
X = sim.data$XMat
Xpred = sim.data$XmatCV
y = sim.data$obsCVBin

d <- rdist.earth(coordinates(obs.grid))
dcross <- rdist.earth(coordinates(obs.grid), coordinates(pred.grid)) # 1600*80
dpred <- rdist.earth(coordinates(pred.grid)) # 80*80

eta.pred <- matrix(NA, nrow = nrow(pred.grid), ncol = niter-burnin) # 80*1600

for(j in 1:ncol(eta.pred)){
  print(j)
  
  # Construct the covariance matrices
  Gamma <- exp(-d/rho.final[j]) 
  Ginv <- solve(Gamma)
  g <- exp(-dcross/rho.final[j])
  Gpred <- exp(-dpred/rho.final[j])
  
  m <- Xpred %*% beta.final[j,] + t(g) %*% Ginv %*% 
    (y - X %*% beta.final[j,])
  V <- s2.final[j] * (Gpred - t(g)%*%Ginv%*%g)
  eta.pred[,j] <- rmvnorm(1, m, V, method = "svd")
}

## Find pointwise posterior means and sds

eta.pred.m <- apply(eta.pred, 1, mean)
eta.pred.sd <- apply(eta.pred, 1, sd)

bin.p = exp(eta.pred.m)/(1+exp(eta.pred.m))

pred.y = sapply(bin.p,rbinom,n=1,size=1)
actual.y = sim.data$obsCVBin

# contingency table
ct = table(pred.y,actual.y)

# accuracy 
sum(diag(ct))/sum(ct)