setwd("C:/Users/user/desktop/jun/SpatioTemp/hw3")
load("munichrents.Rdata")

library(sp)
library(spdep)
library(classInt)
library(fields)
library(MCMCpack)

#-----------------------------------------------------------#
# 1
head(rents)
lm_rents = lm(RentPerM2 ~ . -Location - Room1,data=rents)

# coeff
lm_rents


X = as.matrix(cbind(1,rents[-c(1,3,9)]),nrow=2035)
solve(t(X)%*%X)%*%t(X)%*%as.matrix(rents[c(1)],nrow=2035)
#-----------------------------------------------------------#
# 2 
nb_dist = poly2nb(districts.sp)
coords = coordinates(districts.sp)
png("problem2.png", width = 1000, height = 800)
plot(districts.sp, border = "gray")
plot(nb_dist, coords, pch = 19, cex = 0.6, add = TRUE)
plot(parks.sp,add=TRUE,col="skyblue")
legend("bottomleft", fill = c("skyblue"),
       legend = c("No apartments"),
       bty="n", cex = 1.5, y.intersp = 1.5)
dev.off()

#-----------------------------------------------------------#
# 3
dist_obs = colSums(H)
length(dist_obs)
range(dist_obs)

png("problem3.png", width = 1000, height = 800)

pal <- rev(grey(seq(0,1,length=6))[-1])
q5 <- classIntervals(dist_obs, n = 5, style = "quantile")
col <- findColours(q5, pal)
plot(districts.sp, col = col)
legend("bottomleft", fill = c(attr(col, "palette"),"skyblue"),
       legend = c(names(attr(col, "table")),"No apartments"),
       bty="n", cex = 1.5, y.intersp = 1.5)
plot(parks.sp,add=TRUE,col="skyblue")
dev.off()


#-----------------------------------------------------------#
# 4 

# weight matrix
W = nb2mat(nb_dist,style="B")
Dw = diag(nrow(W))
diag(Dw) = colSums(W)

B = 1e+4 # iterations
p = ncol(X)
m = ncol(H)

# placeholder
beta = matrix(NA, B, p)
eta = matrix(NA, B, m)
sigma2 = matrix(NA, B, 1)
tau2 = matrix(NA, B, 1)


# init
beta[1,] = lm_rents$coeff
eta[1, ] = 1
sigma2[1, ] = 1
tau2[1, ] = 1

# Gibbs sampler
set.seed(1110)

for(i in 2:B){
  start = Sys.time()
  inv_XtX = solve(t(X)%*%X)
  beta[i, ] = rmvnorm(1, mu=inv_XtX%*%t(X)%*%(y-H%*%eta[i-1,]), Sigma=sigma2[i-1,]*inv_XtX)

  inv_HtH = solve(t(H)%*%H/sigma2[i-1,] + (Dw - W)/tau2[i-1,])
  eta[i, ] = rmvnorm(1, mu = inv_HtH%*%t(H)%*%(y-X%*%beta[i,])/sigma2[i-1,], Sigma=inv_HtH)
  eta[i, ] = eta[i,] - mean(eta[i,]) # impose the constraint
  
  sigma2[i, ] = rinvgamma(1,0.001 + n/2, 0.001 + t(y - X%*%beta[i,] - H%*%eta[i,])%*%(y - X%*%beta[i,] - H%*%eta[i,])/2)
  tau2[i, ] = rinvgamma(1,0.001 + (m-1)/2, 0.001 + t(eta[i, ])%*%(Dw- W)%*%eta[i, ]/2)
  if(i ==2) print(paste("ETS: ",round(B*as.vector(Sys.time()-start)/60,3),"mins" ))
  if(i %% 100 == 0){print(i)}
}

save(beta, eta, tau2, sigma2, file = "Gibbs_real_samples.RData")

load("Gibbs_real_samples.RData")

# ACF and trace plot of sigma2 and tau2
png("problem4_tau_ACF.png", width = 1200, height = 600)
par(mfrow=c(1,2))
plot(tau2, type="l", xlab = "iteration", ylab = "tau2", main = "Traceplot of tau2 ")
acf(tau2, xlab = "Lag", ylab = "ACF", main = "ACF plot of tau2 ")
dev.off()

png("problem4_sigma_ACF.png", width = 1200, height = 600)
par(mfrow=c(1,2))
plot(sigma2, type="l", xlab = "iteration", ylab = "sigma2", main = "Traceplot of sigma2 ")
acf(sigma2, xlab = "Lag", ylab = "ACF", main = "ACF plot of sigma2 ")
dev.off()


# posterior mean and C.I

pmeans = matrix(colMeans(beta),nrow=1) # posterior means for beta
colnames(pmeans) = paste0("beta",1:12)
pmeans

CI = apply(beta,2,function(x) quantile(x,probs=c(0.025,0.975)))
colnames(CI) = paste0("beta",1:12)
CI

# map of posterior mean of eta
peta = colMeans(eta)
length(peta)
range(peta)

pal <- rev(grey(seq(0,1,length=6))[-1])
q5 <- classIntervals(peta, n = 5, style = "quantile")
col <- findColours(q5, pal)

png("problem4_post.png", width = 1000, height = 800)

plot(districts.sp, col = col)
legend("bottomleft", fill = c(attr(col, "palette"),"skyblue"),
       legend = c(names(attr(col, "table")),"No apartments"),
       bty="n", cex = 1.5, y.intersp = 1.5)
plot(parks.sp,add=TRUE,col="skyblue")
dev.off()


# map of standard deviation of eta

sd_eta = apply(eta,2,sd)
length(sd_eta)
range(sd_eta)

pal <- rev(grey(seq(0,1,length=6))[-1])
q5 <- classIntervals(sd_eta, n = 5, style = "quantile")
col <- findColours(q5, pal)

png("problem4_sd.png", width = 1000, height = 800)

plot(districts.sp, col = col)
legend("bottomleft", fill = c(attr(col, "palette"),"skyblue"),
       legend = c(names(attr(col, "table")),"No apartments"),
       bty="n", cex = 1.5, y.intersp = 1.5)
plot(parks.sp,add=TRUE,col="skyblue")
dev.off()





