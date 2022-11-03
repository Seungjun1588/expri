########## Ising model example ##########
rm(list=ls())
library(fields)
library(lattice)
library(coda)
library(Rcpp)
library(RcppArmadillo)
library(xtable)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
setwd("C:/Users/JAEWOOPARK/Desktop/Research/Code/Doubly Intractable/review paper code/Ising model example/IsingMCMC")

#####   load saved results   #####
load("Moller.RData")        # call AVM result
load("Murray.RData")        # call Murray's exchange result
load("Liang.RData")         # call Liang DMH result
load("AEX.RData")           # call AEX result
load("Atchade.RData")       # call Atchade's result
load("RRresult.RData")      # call russian roulette result
load("Noisy.RData")         # call Noisy Exchange result
load("NoisyD.RData")        # call Noisy DMH result
load("GoldstandardM.RData") # call gold standard by Murray

#####   Call C++ functions   #####
sourceCpp("IsingMCMC.cpp")   

### Sampling functions ###
# Energy(mat X)
# ProppWilson(int nrow, int ncol, double b)
# Gibb(mat initial, double b, double cycle)
# MH(mat initial, double b, double cycle)
# logZAIS(int noImp, int noTemp, double b, int nrow, int ncol)
# Roulette(int nrow, int ncol, double logZtilde, int NimpSamps, double bprop, double r, double cx, int cmax, int ru, double mu, int NTemp)

### Algorithms ###
# IsingPseudoMCMC(int outer, double initial, double sigma, mat X)
# IsingMoller(int outer, double initial, double sigma, mat X)
# IsingExchange(int outer, double initial, double sigma, mat X)
# IsingAEX(int Niter, int Numaux, double cycle, double t0, int neighbor, vec th, mat thdist, double initial, double sigma, mat X)
# IsingDMH(int outer, double inner, double initial, double sigma, mat X)
# IsingFDMH(int outer, double inner, double initial, double sigma, mat X)
# IsingAtchade(int outer, int inner, double initial, double sigma, vec th, mat X)
# IsingRoulette(int outer, int NimpSamps, int Ntemp, double initial, double sigma, double r, double cx, int cmax, int ru, double mu, mat X)
# IsingNExchange(int outer, int NimpSamps, double initial, double sigma, mat X)
# IsingNDMH(int outer, int inner, int NimpSamps, double initial, double sigma, mat X)

#####   Call batchmean/ESS/Bhattacharyya distance function   #####
source("http://www.stat.psu.edu/~mharan/batchmeans.R")

Bdist = function(a,b,n){
  p = density(a,n=2^15)
  q = density(b,n=2^15)
  sample = c()
  
  for(i in 1:n){
    u = runif(1)
    ind.p = which.min( abs( p$x - u ) )
    ind.q = which.min( abs( q$x - u ) ) 
    sample[i] = sqrt( p$y[ind.p]*q$y[ind.q] )}
  
  result = mean(sample)
  if( result < 1 ){ result = result }else{ 
    result = result - 2*(result - 1) }
  
  return(  -log(result)  )              
}
# sign corrected versio of batch mean
bmsign = function(vals,sign,bs="sqroot",warn=FALSE){
  N <- length(vals)
  if (N<1000)
  {
    if (warn) # if warning
      cat("WARNING: too few samples (less than 1000)\n")
    if (N<10)
      return(NA)
  }
  
  if (bs=="sqroot") 
  {
    b <- floor(sqrt(N)) # batch size
    a <- floor(N/b) # number of batches
  }
  else
    if (bs=="cuberoot") 
    {
      b <- floor(N^(1/3)) # batch size
      a <- floor(N/b) # number of batches
    }
  else # batch size provided
  {
    stopifnot(is.numeric(bs))  
    b <- floor(bs) # batch size
    if (b > 1) # batch size valid
      a <- floor(N/b) # number of batches
    else
      stop("batch size invalid (bs=",bs,")")
  }
  
  Ys <- sapply(1:a,function(k) 
    return(  mean( vals[((k-1)*b+1):(k*b)]*sign[((k-1)*b+1):(k*b)] )/mean(sign[((k-1)*b+1):(k*b)])  ))
  
  muhat <- mean(vals*sign)/mean(sign)
  sigmahatsq <- b*sum((Ys-muhat)^2)/(a-1)
  
  bmse <- sqrt(sigmahatsq/N)
  
  return(list(est=muhat,se=bmse))
}

##### run step #####

### Simulate the data ### 7
set.seed(7)
N <- 10                     # number of lattice size
X <- ProppWilson(N,N,0.2)   # assume it is true data 
n <- 20000                  # number of MCMC RUN

#### Conduct Moller's Auxiliary method ###
## log likelihood for MPLE ##
pseudoLik <- function(parameter){
  beta <- parameter
  work <- cbind(0, rbind(0, X, 0), 0); c <- 0
  for(i in 2:(N+1)){
    for(j in 2:(N+1)){
      p <- exp(2*beta*(work[i,j-1]+work[i,j+1]+work[i-1,j]+work[i+1,j]))
      p <- p/(1+p)
      c <- c + dbinom((work[i,j]+1)/2,1,p,log=TRUE)}}
  return(c)
}

z <- optim(0.1,pseudoLik,control=list(fnscale=-1), method = "BFGS")
MPLE <- z$par; MPLE

ptm <- proc.time()
Moller <- IsingMoller(n,MPLE,0.1,X)
ptme <- proc.time() 

par(mfrow=c(1,2))
hist(Moller)
ts.plot(Moller)

bm(Moller)
HPDinterval(as.mcmc(Moller), prob = 0.95)
length(unique(Moller))/n
ess(Moller)
Mollersummary <- c(bm(Moller)$est, bm(Moller)$se,HPDinterval(as.mcmc(Moller), prob = 0.95)[1],
                   HPDinterval(as.mcmc(Moller), prob = 0.95)[2],ess(Moller),length(unique(Moller))/n,(ptme-ptm)[1])

save(Mollersummary,Moller,file="Moller.RData")

#### Conduct Murray's Exchange ###

ptm <- proc.time()
Murray <- IsingExchange(n,MPLE,0.1,X)
ptme <- proc.time() 

par(mfrow=c(1,2))
hist(Murray)
ts.plot(Murray)

bm(Murray)
HPDinterval(as.mcmc(Murray), prob = 0.95)
length(unique(Murray))/n
ess(Murray)
Murraysummary <- c(bm(Murray)$est, bm(Murray)$se,HPDinterval(as.mcmc(Murray), prob = 0.95)[1],
                   HPDinterval(as.mcmc(Murray), prob = 0.95)[2],ess(Murray),length(unique(Murray))/n,(ptme-ptm)[1])

save(Murraysummary,Murray,file="Murray.RData")

#### Conduct Liang's DMH ###

ptm <- proc.time()
Liang <- IsingDMH(n,10,MPLE,0.1,X)
ptme <- proc.time() 

par(mfrow=c(1,2))
hist(Liang)
ts.plot(Liang)

bm(Liang)
HPDinterval(as.mcmc(Liang), prob = 0.95)
length(unique(Liang))/n
ess(Liang)
Liangsummary <- c(bm(Liang)$est, bm(Liang)$se,HPDinterval(as.mcmc(Liang), prob = 0.95)[1],
                  HPDinterval(as.mcmc(Liang), prob = 0.95)[2],ess(Liang),length(unique(Liang))/n,(ptme-ptm)[1])

save(Liangsummary, Liang,file="Liang.RData")


#### Conduct AEX ###
load("FLiang.RData")

# step 1. Conduct Liang's Fractional DMH to get auxiliary parameters # 
m <- 100   # number of auxiliary parameters
aux.par <- rep(0,m)
N1 <- 5500

ptm1 <- proc.time()
FLiang <- IsingFDMH(N1,10,MPLE,0.1,X)  # multiply 0.5 in c code
ptme1 <- proc.time() 

par(mfrow=c(1,2))
hist(FLiang)
ts.plot(FLiang)
save(FLiang,file="FLiang.RData")

FLiang <- FLiang[501:5500]                               # burn in 500
stand <-  (FLiang-min(FLiang))/(max(FLiang)-min(FLiang)) # standardized  
stand <- unique(stand)                                   # only take unique components
dmat <- rdist(stand)                                     # distance mat

# choose auxiliary parameters through min max procedure
ind = 1; A = 1; Ac = 2:length(stand)
aux.par[1] = stand[ind]

ind = which.max( dmat[,A] )
A = c(A,ind)
Ac = Ac[-which(Ac==ind)]
aux.par[2] = stand[ind]


for(i in 3:m){
  dummy = max( apply( dmat[,A] , 1, min )[Ac] )
  ind = which( dmat[,A] == dummy  )
  if(ind < dim(dmat)[1]){ ind = ind }else{ ind = ind-floor( ind/dim(dmat)[1] )*dim(dmat)[1] }
  A = c(A,ind)
  Ac = Ac[-which(Ac==ind)]
  aux.par[i] = stand[ind]
}

dist.aux.par = rdist(aux.par)  # distance matrix for aux.par (for standardized version)
aux.par = (max(FLiang)-min(FLiang))*aux.par + min(FLiang)  


# step 2. Run AEX # 

Numaux = 20000
t0 = 20000
neighbor = 30
sigma = 0.1

ptm <- proc.time()
res = IsingAEX(n, Numaux, 1, t0, neighbor, aux.par, dist.aux.par, MPLE,sigma, X)
ptme <- proc.time() 
AEX = res$par

par(mfrow=c(1,2))
hist(AEX)
ts.plot(AEX)

bm(AEX)
HPDinterval(as.mcmc(AEX), prob = 0.95)
length(unique(AEX))/n
ess(AEX)
AEXsummary <- c(bm(AEX)$est, bm(AEX)$se,HPDinterval(as.mcmc(AEX), prob = 0.95)[1],
                HPDinterval(as.mcmc(AEX), prob = 0.95)[2],ess(AEX),length(unique(AEX))/n,(ptme1-ptm1)[1]+(ptme-ptm)[1])

save(AEXsummary, AEX,file="AEX.RData")

#### Conduct Atchade's Adaptive method ###
th <- runif(100,0,1)
# for(i in 1:100){th[i] = Liang[200*i]} # possible add (ptme1-ptm1)[1]

ptm <- proc.time()
Atchade <- IsingAtchade(n,1,MPLE,0.1,th,X)
ptme <- proc.time() 

par(mfrow=c(1,2))
hist(Atchade)
ts.plot(Atchade)

bm(Atchade)
HPDinterval(as.mcmc(Atchade), prob = 0.95)
length(unique(Atchade))/n
ess(Atchade)
Atchadesummary <- c(bm(Atchade)$est, bm(Atchade)$se,HPDinterval(as.mcmc(Atchade), prob = 0.95)[1],
                    HPDinterval(as.mcmc(Atchade), prob = 0.95)[2],ess(Atchade),length(unique(Atchade))/n,(ptme-ptm)[1])

save(Atchadesummary, Atchade,file="Atchade.RData")

#### Conduct Russian Roulette (pseudomarginal) ###
## 4 and half hour ## 
Nsamps = n         # no of iterations in chain
NimpSamps = 1000   # no importance samples
cmax = 50          # maximum no of terms in series before resorting to poisson truncation
cx = 0.4           # kappa = 1 - cx*(Z_est/Z_tilde)
r = 0.6            # roulette parameter
russ = 1           # set to 1 to use russian roulette, 0 to use poisson truncation
poiss_mu = 1       # mean of poisson distribution
stepb = 0.1        # sd of proposal
noTemp = 2000      # number of temperature levels in AIS

ptm <- proc.time()
out <- IsingRoulette(n, NimpSamps, noTemp, MPLE,stepb, r, cx, cmax, russ, poiss_mu, X)
ptme <- proc.time() 

RussianR<- out[1,] 
si <- out[2,]

par(mfrow=c(1,2))
hist(RussianR)
ts.plot(RussianR)

# ean and std of samples and absolute corrected version
ExpS_noCorr = mean(RussianR)
ExpS_Corr = mean(RussianR*si)/mean(si)
StdS_noCorr = sd(RussianR)
Std_Corr = sqrt(mean(RussianR^2*si)/mean(si) - ExpS_Corr^2)

bmsign(RussianR,si)
HPDinterval(as.mcmc(RussianR), prob = 0.95)
length(unique(RussianR))/n
ess(RussianR)
RussianRsummary <- c(bmsign(RussianR,si)$est, bmsign(RussianR,si)$se,HPDinterval(as.mcmc(RussianR), prob = 0.95)[1],
                     HPDinterval(as.mcmc(RussianR), prob = 0.95)[2],ess(RussianR),length(unique(RussianR))/n,((ptme-ptm)[1]/8))

save(RussianRsummary,out,file="RRresult.RData")

#### Conduct Noisy exchange ###
# takes around 12min #
ptm <- proc.time()
Noisy <- IsingNExchange(n, 100, MPLE,0.1, X)
ptme <- proc.time() 

par(mfrow=c(1,2))
hist(Noisy)
ts.plot(Noisy)

bm(Noisy)
HPDinterval(as.mcmc(Noisy), prob = 0.95)
length(unique(Noisy))/n
ess(Noisy)
Noisysummary <- c(bm(Noisy)$est, bm(Noisy)$se,HPDinterval(as.mcmc(Noisy), prob = 0.95)[1],
                  HPDinterval(as.mcmc(Noisy), prob = 0.95)[2],ess(Noisy),length(unique(Noisy))/n,((ptme-ptm)[1]/8))

save(Noisysummary,Noisy,file="Noisy.RData")

#### Conduct Noisy DMH ###
# takes around 12min #
ptm <- proc.time()
NoisyD <- IsingNDMH(n,10,100,MPLE,0.1,X)
ptme <- proc.time() 

par(mfrow=c(1,2))
hist(NoisyD)
ts.plot(NoisyD)

bm(NoisyD)
HPDinterval(as.mcmc(NoisyD), prob = 0.95)
length(unique(NoisyD))/n
ess(NoisyD)
NoisyDsummary <- c(bm(NoisyD)$est, bm(NoisyD)$se,HPDinterval(as.mcmc(NoisyD), prob = 0.95)[1],
                   HPDinterval(as.mcmc(NoisyD), prob = 0.95)[2],ess(NoisyD),length(unique(NoisyD))/n,((ptme-ptm)[1]/8))

save(NoisyDsummary,NoisyD,file="NoisyD.RData")


#### Goldstandard through Murray's Exchange ###

ptm <- proc.time()
GoldstandardM <- IsingExchange(1010000,MPLE,0.1,X)
proc.time() - ptm
save(GoldstandardM,file="GoldstandardM.RData")


# burn in 10000 and choose samples 
GoldstandardM <- GoldstandardM[10001:1010000] 
Gold <- c()
for(i in 1:10000){ Gold[i] = GoldstandardM[100*i] }

par(mfrow=c(1,2))
hist(Gold)
ts.plot(Gold)

bm(Gold)
HPDinterval(as.mcmc(Gold), prob = 0.95)
ess(Gold)
Goldsummary <- c(bm(Gold)$est, bm(Gold)$se,HPDinterval(as.mcmc(Gold), prob = 0.95)[1],
                 HPDinterval(as.mcmc(Gold), prob = 0.95)[2],ess(Gold),NA,NA)


save.image("100lattice020.RData")
#### Comparison analysis ###

# Comparison with Gold standard #
# In exact MCMC (Liang, NoisyDMH) are tuned epsilon < 0.01 #
Bdist(Murray,Gold,100000)
Bdist(Moller,Gold,100000)
Bdist(Atchade,Gold,100000)
Bdist(RussianR,Gold,100000)
Bdist(Liang,Gold,100000)
Bdist(NoisyD,Gold,100000)

pdf("Isingdensity.pdf")
plot(density(Gold),main=expression(theta==0.2),xlab= "parameter space",ylim=c(0,6.2),lwd=3,col="goldenrod3")
lines(density(Moller),lwd=3,col=1)
lines(density(Murray),lwd=3,col=2)
lines(density(Liang),lwd=3,col=3)
lines(density(AEX),lwd=3,col=4)
lines(density(Atchade),lwd=3,col=5)
lines(density(RussianR),lwd=3,col=6)
lines(density(NoisyD),lwd=3,col=7)

legend("topright",c("Gold","AVM","Exchange","DMH","AEX","ALR","RussianR","NoisyDMH"),
       col=c("goldenrod3",1,2,3,4,5,6,7),lwd=rep(3,8))
dev.off()


pdf("Isingtraceplot.pdf")
par(mfrow=c(2,3))
ts.plot(Moller,ylab="AVM")
ts.plot(Murray,ylab="Exchange")
ts.plot(Liang,ylab="DMH")
ts.plot(AEX)
ts.plot(Atchade)
ts.plot(RussianR)
ts.plot(NoisyD,ylab="NoisyDMH")
dev.off()

pdf("Isingacf.pdf")
par(mfrow=c(2,4))
acf(Moller,main="AVM")
acf(Murray,main="Exchange")
acf(Liang,main="DMH")
acf(AEX,main="AEX")
acf(Atchade,main="Atchade")
acf(RussianR,main="RussianR")
acf(NoisyD,main="NoisyDMH")
dev.off()

table1 = rbind(Mollersummary, Murraysummary, Liangsummary, AEXsummary,
               Atchadesummary, RussianRsummary, NoisyDsummary,Goldsummary )

rownames(table1)=c("AVM","Exchange","DMH","AEX","Atchade","RussianR","NoisyDMH","Gold")
colnames(table1)=c("Mean","MCSE","95%HPD","95%HPD","ESS","Acc","Time(second)")
table1 = round(table1,2)
table1
xtable(table1,auto=TRUE)



ESSperT = function(x){round(x[5]/x[7],3)}
ESSperT(Mollersummary)
ESSperT(Murraysummary)
ESSperT(Liangsummary)
ESSperT(AEXsummary)
ESSperT(Atchadesummary)
ESSperT(RussianRsummary)
ESSperT(NoisyDsummary)



#### Comparison analysis for lattice 045 ###

load("100lattice045.RData")
RussianR<- out[1,] 
si <- out[2,]

Gold <- c()
for(i in 1:10000){ Gold[i] = GoldstandardM[100*i] }
par(mfrow=c(1,2))
hist(Gold)
ts.plot(Gold)

bm(Gold)
HPDinterval(as.mcmc(Gold), prob = 0.95)
ess(Gold)
Goldsummary <- c(bm(Gold)$est, bm(Gold)$se,HPDinterval(as.mcmc(Gold), prob = 0.95)[1],
                 HPDinterval(as.mcmc(Gold), prob = 0.95)[2],ess(Gold),NA,NA)



pdf("Isingdensity043.pdf")
plot(density(Gold),main=expression(theta==0.43),xlab= "parameter space",ylim=c(0,8.5),lwd=3,col="goldenrod3")
lines(density(Moller),lwd=3,col=1)
lines(density(Murray),lwd=3,col=2)
lines(density(Liang),lwd=3,col=3)
lines(density(AEX),lwd=3,col=4)
lines(density(Atchade),lwd=3,col=5)
lines(density(RussianR),lwd=3,col=6)
lines(density(NoisyD),lwd=3,col=7)
legend("topright",c("Gold","AVM","Exchange","DMH","AEX","ALR","RussianR","NoisyDMH"),
       col=c("goldenrod3",1,2,3,4,5,6,7),lwd=rep(3,8))
dev.off()

table1 = rbind(Mollersummary, Murraysummary, Liangsummary, AEXsummary,
               Atchadesummary, RussianRsummary, NoisyDsummary,Goldsummary )

rownames(table1)=c("AVM","Exchange","DMH","AEX","Atchade","RussianR","NoisyDMH","Gold")
colnames(table1)=c("Mean","MCSE","95%HPD","95%HPD","ESS","Acc","Time(second)")
table1 = round(table1,2)
table1
xtable(table1,auto=TRUE)



ESSperT = function(x){round(x[5]/x[7],3)}
ESSperT(Mollersummary)
ESSperT(Murraysummary)
ESSperT(Liangsummary)
ESSperT(AEXsummary)
ESSperT(Atchadesummary)
ESSperT(RussianRsummary)
ESSperT(NoisyDsummary)



##############################################################################################################
## FINAL PLOT
##############################################################################################################
# density plot
pdf("isingdensity.pdf",width=8,height=4)
par(mfrow=c(1,2), mar=c(4,4,2,2)) 
load("100lattice020.RData")
plot(density(Gold),main=expression(theta==0.2),xlab= "parameter space",ylim=c(0,6.2),lwd=3,col="goldenrod3")
lines(density(Moller),lwd=3,col=1)
lines(density(Murray),lwd=3,col=2)
lines(density(Liang),lwd=3,col=3)
lines(density(AEX),lwd=3,col=4)
lines(density(Atchade),lwd=3,col=5)
lines(density(RussianR),lwd=3,col=6)
lines(density(NoisyD),lwd=3,col=7)
legend("topright",c("Gold","AVM","Exchange","DMH","AEX","ALR","RussianR","NoisyDMH"),
       col=c("goldenrod3",1,2,3,4,5,6,7),lwd=rep(3,8))


load("100lattice045.RData")
RussianR<- out[1,] 
si <- out[2,]
Gold <- c()
for(i in 1:10000){ Gold[i] = GoldstandardM[100*i] }
bm(Gold)
HPDinterval(as.mcmc(Gold), prob = 0.95)
ess(Gold)
Goldsummary <- c(bm(Gold)$est, bm(Gold)$se,HPDinterval(as.mcmc(Gold), prob = 0.95)[1],
                 HPDinterval(as.mcmc(Gold), prob = 0.95)[2],ess(Gold),NA,NA)

plot(density(Gold),main=expression(theta==0.43),xlab= "parameter space",ylim=c(0,8.5),lwd=3,col="goldenrod3")
lines(density(Moller),lwd=3,col=1)
lines(density(Murray),lwd=3,col=2)
lines(density(Liang),lwd=3,col=3)
lines(density(AEX),lwd=3,col=4)
lines(density(Atchade),lwd=3,col=5)
lines(density(RussianR),lwd=3,col=6)
lines(density(NoisyD),lwd=3,col=7)
legend("topright",c("Gold","AVM","Exchange","DMH","AEX","ALR","RussianR","NoisyDMH"),
       col=c("goldenrod3",1,2,3,4,5,6,7),lwd=rep(3,8))

dev.off()


# boxplot
pdf("isingdensity.pdf")
par(mfrow=c(2,1), mar=c(2,2,2,2))
load("100lattice020.RData")
boxplot(Gold,Moller,Murray,Liang,AEX,Atchade,RussianR,NoisyD,outline=FALSE,
        names=c("Gold","AVM","Exchange","DMH","AEX","ALR","Russian","NoisyD"),main=expression(theta==0.2))
abline(h=0.2,col=2,lty=3,lwd=3)


load("100lattice045.RData")
RussianR<- out[1,] 
si <- out[2,]
Gold <- c()
for(i in 1:10000){ Gold[i] = GoldstandardM[100*i] }
bm(Gold)
HPDinterval(as.mcmc(Gold), prob = 0.95)
ess(Gold)
Goldsummary <- c(bm(Gold)$est, bm(Gold)$se,HPDinterval(as.mcmc(Gold), prob = 0.95)[1],
                 HPDinterval(as.mcmc(Gold), prob = 0.95)[2],ess(Gold),NA,NA)
boxplot(Gold,Moller,Murray,Liang,AEX,Atchade,RussianR,NoisyD,outline=FALSE,
        names=c("Gold","AVM","Exchange","DMH","AEX","ALR","Russian","NoisyD"),main=expression(theta==0.43))
abline(h=0.43,col=2,lty=3,lwd=3)
dev.off()



