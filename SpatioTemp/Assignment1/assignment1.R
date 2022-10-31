rm(list=ls())

setwd("C:/Users/user/Desktop/jun/SpatioTemp")


#-------------------------------------------------------------#
# 1-(a)

# 1-(b)

samp_mat = matrix(c(1,0.5,0.5,1),2)

sigma = samp_mat
sim_y = function(mu,sigma){
  L =t(chol(sigma)) # lower triangle matrix
  n = length(mu)
  Z = rnorm(n,0,1)
  Y = mu + L%*%Z
  return(Y)
}


#1-(c)

x = seq(0, 1, length = 500) # fine grid
d = as.matrix(dist(x)) # create distance matrix



# We consider an exponential cov function.
par(mfrow=c(1,3))

sigma1 <- 10*exp(-d/0.5) # cov matrix
samp1 = sim_y(mu=rep(0,500),sigma1)
samp2 = sim_y(mu=rep(0,500),sigma1)
samp3 = sim_y(mu=rep(0,500),sigma1)
samples1= cbind(samp1,samp2,samp3)
matplot(x, samples1, type = "l",main = expression(rho*"=0.5"))


sigma2 <- 10*exp(-d/10) # cov matrix
samp1 = sim_y(mu=rep(0,500),sigma2)
samp2 = sim_y(mu=rep(0,500),sigma2)
samp3 = sim_y(mu=rep(0,500),sigma2)
samples2= cbind(samp1,samp2,samp3)
matplot(x, samples2, type = "l",main = expression(rho*"=10"))

samp1 = sim_y(mu=rep(0,500),sigma2)
samp2 = sim_y(mu=rep(3,500),sigma2)
samp3 = sim_y(mu=rep(5,500),sigma2)
samples2= cbind(samp1,samp2,samp3)
matplot(x, samples2, type = "l",main = expression(mu*"=0"*", "*mu*"=3"*", "*mu*"=5"))


#-------------------------------------------------------------#
# 2
library(sp)
library(gstat)
library(fields)
library(classInt)
library(maps)

# 2-(a)
load("hw1/CAtemps.RData")


head(CAtemp)

data = data.frame(cbind(coordinates(CAtemp),CAtemp$avgtemp,CAtemp$elevation))
sum(is.na(data)) # no missing data

colnames(data) = c("lon","lat","avgtemp","elevation")
fit = lm(avgtemp ~ lon + lat + elevation, data=data)
summary(fit)
coeff = fit$coefficients  # answer 
coeff
res = data$avgtemp - fit$fitted.values

# plot.point.ref(CAtemp, res)


ploteqc <- function(spobj, z, breaks, ...){
  pal <- tim.colors(length(breaks)-1)
  fb <- classIntervals(z, n = length(pal), 
                       style = "fixed", fixedBreaks = breaks)
  col <- findColours(fb, pal)
  plot(spobj, col = col, ...)
  image.plot(legend.only = TRUE, zlim = range(breaks), col = pal)
}

range(res)
breaks = seq(-7, 7, length.out = 10)
ploteqc(CAtemp, res, breaks, pch = 19) # answer
map("county", region = "california", add = TRUE)
title(main = "Average Annual Temperatures, 1961-1990, Degrees F")


# 2-(b)

CAtemp$res = res

vg = variogram(res ~ 1, data = CAtemp, width=40)
plot(vg, xlab = "Distance", ylab = "Semi-variogram estimate", width=5) # What happens when we 
# change the bin width?

# vgangle = variogram(res ~ 1, data = CAtemp, alpha = c(0, 45, 90, 135))
# plot(vgangle, xlab = "Distance", ylab = "Semi-variogram estimate")
# # I think there are no big differences..

vgmodel = vgm(1, "Exp", 200, 0.05)
fitvg = fit.variogram(vg, vgmodel)
plot(vg, fitvg, xlab = "Distance", ylab = "Semi-variogram estimate")


print(fitvg)
s2.hat = fitvg$psill[2]
rho.hat = fitvg$range[2]
tau2.hat = fitvg$psill[1]
c(s2.hat,rho.hat,tau2.hat)


# 2-(c)
head(data)
data2 = cbind(coordinates(CAtemp),CAtemp$elevation) 
colnames(data2) = c("lon","lat","elevation")
X = as.matrix(cbind(1,data2))

dist = rdist.earth(coordinates(CAtemp))
covmat = exp(-dist/rho.hat)*s2.hat
diag(covmat) = s2.hat
nugget = diag(nrow(dist))*tau2.hat
Sigma = covmat + nugget
Sigma_inv = solve(covmat + nugget)
y = CAtemp$avgtemp
beta_gls = solve(t(X)%*% Sigma_inv %*% X) %*% t(X) %*% Sigma_inv %*% y

beta_gls

#2-(d)
dd = rdist.earth(coordinates(CAtemp),coordinates(CAgrid))
gamma = s2.hat * exp(-dd/rho.hat) 
Xs = cbind(1,coordinates(CAgrid),CAgrid$elevation)

ypred = Xs %*% beta_gls + t(gamma) %*% solve(Sigma) %*%(y - X %*% beta_gls)

range(ypred)
breaks = seq(30, 75, length.out = 10)
ploteqc(CAgrid, ypred, breaks, pch = 19) # answer
map("county", region = "california", add = TRUE)
points(coordinates(CAtemp),pch=4,cex=1)
title(main = "Mean estimate Annual Temperatures, 1961-1990, Degrees F")


b = t(Xs) - t(X) %*% solve(Sigma) %*% gamma
vpred = s2.hat - diag(t(gamma) %*% Sigma_inv%*%gamma + t(b) %*% solve(t(X) %*% Sigma_inv %*% X ) %*% b)


sepred = sqrt(vpred)

range(sepred)
breaks = seq(range(sepred)[1],range(sepred)[2], length.out = 10)
ploteqc(CAgrid, sepred, breaks, pch = 19) # answer
map("county", region = "california", add = TRUE)
points(coordinates(CAtemp),pch=4,cex=1)
title(main = "Standard error Annual Temperatures, 1961-1990, Degrees F")




