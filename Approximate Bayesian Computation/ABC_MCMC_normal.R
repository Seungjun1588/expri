###
### mean 5, sd =known normal simulation
###

#-----------------------------------------------------------#
# data 
rm(list=ls())
set.seed(123450)
n=100
true_mean = 30 # unknown
sd = 12  # known
y = rnorm(n,true_mean,sd)
mean(y)
#-----------------------------------------------------------#
# prior 
theta0 = 20 # prior
sigma0 = 10 # prior

#-----------------------------------------------------------#
distance1 = function(x,y){
  return((sum(x)-sum(y))^2)
}

distance2 = function(x,y){
  return(abs(mean(x)-mean(y)))
}

distance3 = function(x,y){
  return((mean(x)-mean(y))^2)
}
#-----------------------------------------------------------#
# hyper pararmeter


#-----------------------------------------------------------#
# MCMC



MCMC = function(N,step_size){
  result1 =c()
  
  # kernel of likelihood
  likeli = function(y,mu,sigma){
    summ = sum((y-mu)^2)
    return(exp((-1/(2*sigma^2))*summ))
  }
  
  # kernel of prior
  prior = function(new_theta,theta0,tau0){
    return(exp((-1/(2*tau0^2))*(new_theta-theta0)))
  }
  
  theta = 20 # initial value
  
  for(i in 1:N){
    new_theta = rnorm(1,theta,step_size)
    ratio = (likeli(y=y,mu=new_theta,sigma=sd)*prior(new_theta=new_theta,theta0=theta0,tau0=sigma0))/(likeli(y=y,mu=theta,sigma=sd)*prior(new_theta=theta,theta=theta0,tau0=sigma0))
    
    
    if( (runif(1) < ratio )){
      result1[i] = new_theta
      theta = new_theta
    }
    else{
      result1[i] = theta
    }
  }
  
  burnin = 0.2*N
  result1 = result1[-(1:burnin)]
  return(result1)
}
result1 = MCMC(N=1e+06,step_size=1)


ts.plot(result1,main="MCMC method")
acf(result1,main="MCMC method")
length(result1)
mean(result1)
sd(result1)  

#-----------------------------------------------------------#
# ABC(mean)

ABC = function(N,epsilon){
  result2 = c()
  
  for(i in 1:N){
    theta = rnorm(1,mean=theta0,sd=sigma0) # prior
    z = rnorm(n,theta,sd) # sd는 이미 알고 있는 상황
    
    if( distance2(y,z) < epsilon){
      result2= c(result2,theta)
    }
  }
  return(result2)
}

result2 = ABC(N=1e+06,epsilon=0.1)
length(result2)

ts.plot(result2,main="ABC method")
mean(result2)
sd(result2)

#-----------------------------------------------------------#
# MCMC_ABC
MCMC_ABC = function(N,step_size,epsilon){
  result3 =c()
  
  # kernel of likelihood
  likeli = function(y,mu,sigma){
    summ = sum((y-mu)^2)
    return(exp((-1/(2*sigma^2))*summ))
  }
  
  # kernel of prior
  prior = function(new_theta,theta0,tau0){
    return(exp((-1/(2*tau0^2))*(new_theta-theta0)))
  }
  
  theta = result2[length(result2)] 
  
  for(i in 1:N){
    new_theta = rnorm(1,theta,step_size)
    z = rnorm(n,new_theta,sd) # sd는 이미 안다고 가정
    ratio = prior(new_theta=new_theta,theta0=theta0,tau0=sigma0)/prior(new_theta=theta,theta0=theta0,tau0=sigma0)
    
    
    if( (runif(1) < ratio) & (distance2(y,z) < epsilon)){
      result3[i] = new_theta
      theta = new_theta
    }
    else{
      result3[i] = theta
    }
  }
  burnin = 0.2*N
  result3 = result3[-(1:burnin)]
  return(result3)
}

result3 = MCMC_ABC(N=1e+06,step_size=5,epsilon=0.1)

ts.plot(result3,main="MCMC_ABC method")
acf(result3,main="MCMC_ABC method")

length(result3)
mean(result3)
sd(result3)
#-----------------------------------------------------------#

# posterior dist
# conjugacy를 만족하기 때문에 정확한 posterior 구할 수 있다.

n=length(y)
N = 1e+05
mu_n = ((1/sigma0^2)/(1/sigma0^2+n/sd^2))*theta0 + ((n/sd^2)/(1/sigma0^2 + n/sd^2))*mean(y)
tau_sq = 1/(1/sigma0^2 +n/sd^2)
c("post_mean"=mu_n,"post_var"=tau_sq)
post_sample = rnorm(N,mu_n,tau_sq)
mean(post_sample)
sd(post_sample)
plot(density(post_sample),lwd=2,col="black",main="Density estimation",xilm=c(25,35),ylim=c(0,0.4))


# density curve

lines(density(result1),col="blue",lwd=2)
lines(density(result2),col="red",lwd=2)
lines(density(result3),col="green",lwd=2)
legend("topleft",legend=c("True posterior","MCMC","ABC","MCMC_ABC"),fill=c("black","blue","red","green"))

length(result1)
length(result2)
length(result3)
#--------------------------------------------------------------#
# ABC epsilon 다르게 해보기
niter = 1e+05
ABC_0.05 = ABC(N=niter,epsilon=0.05)
ABC_0.1 = ABC(N=niter,epsilon=0.1)
ABC_0.5 = ABC(N=niter,epsilon=0.5)
ABC_1 = ABC(N=niter,epsilon=1)

length(ABC_0.05)
length(ABC_1)

plot(density(post_sample),lwd=2,col="black",main="ABC Density estimation with diff. epsilon",xilm=c(25,35),ylim=c(0,0.4))
lines(density(ABC_0.05),col="blue",lwd=2)
lines(density(ABC_0.1),col="red",lwd=2)
lines(density(ABC_0.5),col="green",lwd=2)
lines(density(ABC_1),col="yellow",lwd=2)
lines(density(post_sample),col="black",lwd=2)
legend("topleft",legend=c("True posterior","0.05","0.1","0.5","1"),fill=c("black","blue","red","green","yellow"))
#--------------------------------------------------------------#
# ABC iter 다르게 해보기


niter = 1e+04
ABC_e4 = ABC(N=niter,epsilon=0.05)

niter = 1e+05
ABC_e5 = ABC(N=niter,epsilon=0.05)
niter = 1e+05

niter = 1e+06
ABC_e6 = ABC(N=niter,epsilon=0.1)

niter = 1e+07
ABC_e7 = ABC(N=niter,epsilon=0.5)


plot(density(post_sample),lwd=2,col="black",main="ABC Density estimation with diff. epsilon",xilm=c(25,35),ylim=c(0,0.4))
lines(density(ABC_0.05),col="blue",lwd=2)
lines(density(ABC_0.1),col="red",lwd=2)
lines(density(ABC_0.5),col="green",lwd=2)
lines(density(ABC_1),col="yellow",lwd=2)
lines(density(post_sample),col="black",lwd=2)
legend("topleft",legend=c("True posterior","0.05","0.1","0.5","1"),fill=c("black","blue","red","green","yellow"))


#--------------------------------------------------------------#
# different step size
niter = 1e+06
MCMC1 = MCMC(N=niter,step_size=1)
MCMC5 = MCMC(N=niter,step_size=5)
MCMC10 = MCMC(N=niter,step_size=10)


plot(density(post_sample),lwd=2,col="black",main="MCMC Density estimation with diff. step size",xilm=c(25,35),ylim=c(0,0.4))
lines(density(MCMC1),col="blue",lwd=2)
lines(density(MCMC5),col="red",lwd=2)
lines(density(MCMC10),col="green",lwd=2)
lines(density(post_sample),col="black",lwd=2)
legend("topleft",legend=c("True posterior","1","5","10"),fill=c("black","blue","red","green"))
#--------------------------------------------------------------#
# different data size

niter = 1e+05
n=100
MCMC100 = MCMC(N=niter,step_size=1)
n=500
MCMC500 = MCMC(N=niter,step_size=1)
n=10000
MCMC10000 = MCMC(N=niter,step_size=1)


plot(density(post_sample),lwd=2,col="black",main="MCMC Density estimation with diff. data size",xilm=c(25,35),ylim=c(0,0.4))
lines(density(MCMC1),col="blue",lwd=2)
lines(density(MCMC5),col="red",lwd=2)
lines(density(MCMC10),col="green",lwd=2)
lines(density(post_sample),col="black",lwd=2)
legend("topleft",legend=c("True posterior","1","5","10"),fill=c("black","blue","red","green"))

