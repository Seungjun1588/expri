###
### p =0.5 binomial simulation
###
#-----------------------------------------------------------#
# data ~ B(p=?)
#set.seed(11211)
rm(list=ls())
n =100
y = rbinom(n,1,p=0.5)
mean(y)

s = sum(y) # success
f = sum(y) # failure
#-----------------------------------------------------------#
# prior ~ Beta(1,1)
a0 = 1 # prior
b0 = 1 # prior

#-----------------------------------------------------------#
distance1 = function(x,y){
  return(sum((x-y)^2))
}

distance2 = function(x,y){
  return(abs(mean(x)-mean(y)))
}

distance3 = function(x,y){
  return(abs(mean(x)-mean(y)))
}
#-----------------------------------------------------------#
# hyper pararmeter
N = 1e+05
epsilon = 0.1
step_size = 0.5 # proposal 


#-----------------------------------------------------------#
# MCMC

result1 =c()


theta = 0.3 # intial value
for(i in 1:N){
  new_theta = runif(1)
  ratio = ((new_theta^(s+a0-1))*(1-new_theta^(f+b0-1)))/((theta^(s+a0-1))*(1-theta^(f+b0-1)))

  if( (0 < new_theta ) &&(new_theta < 1)&&(runif(1) < ratio )){
    result1[i] = new_theta
    theta = new_theta
  }
  else{
    result1[i] = theta
  }
}

burnin = 0.2*N
result1 = result1[-(1:burnin)]

ts.plot(result1,main="MCMC method")
acf(result1,main="MCMC method")
length(result1)
mean(result1)
sd(result1)

#-----------------------------------------------------------#
# ABC(mean)
result2 = c()
epsilon=0.5
for(i in 1:N){
  new_theta = runif(1)
  z = rbinom(n,1,new_theta)
  
  if( (distance2(z,y) < epsilon )){
    result2[i] = new_theta
    theta = new_theta
  }
  else{
    result2[i] = theta
  }
}
length(result2)
acf(result2)

ts.plot(result2,main="ABC method")
acf(result2,main="ABC method")
mean(result2)
sd(result2)

#-----------------------------------------------------------#
# MCMC_ABC
result3 =c()

theta = result2[length(result2)] 

for(i in 1:N){
  new_theta = rbeta(1,1,1)
  z = rbinom(n,1,new_theta)
  ratio = ((new_theta^(s))*(1-new_theta^(f)))/((theta^(s))*(1-theta^(f)))
  
  
  if( (runif(1) < ratio ) && (distance2(z,y) < epsilon )){
    result3[i] = new_theta
    theta = new_theta
  }
  else{
    result3[i] = theta
  }
}

ts.plot(result3,main="MCMC_ABC method")
acf(result3,main="MCMC_ABC method")
burnin = 0.2*N
result3 = result3[-(1:burnin)]


length(result3)
mean(result3)
sd(result3)
#-----------------------------------------------------------#

# posterior dist
# conjugacy를 만족하기 때문에 정확한 posterior 구할 수 있다.

n=length(y)
a_n = s+a0
b_n = f + b0
c("post_alpha"=a_n,"post_beta"=b_n)
post_sample = rbeta(N,a_n,b_n)
mean(post_sample)
sd(post_sample)
plot(density(post_sample),lwd=2,col="black",main="Density estimation")


# density curve

lines(density(result1),col="blue",lwd=2)
lines(density(result2),col="red",lwd=2)
lines(density(result3),col="green",lwd=2)
legend("topright",legend=c("Posterior","MCMC","ABC","MCMC_ABC"),fill=c("black","blue","red","green"))

## ABC가 계속 편향이 있는 것처럼 나오는 이유
# ABC는 prior가 unif이기 때문에 나오는 차이이다. 

