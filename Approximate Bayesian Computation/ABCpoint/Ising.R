#Proppwilson 을 보라고 하신 이유는 이해를 못했음...


rm(list=ls())
setwd("C:/Users/user/Desktop/승준/연구/Sx_DNN")

#-------------------------------------------------------------------------------#

ising.mh=function(nx = 50, ny=20, beta = .2, istop = 100){
  x <- matrix(0,nrow=nx+2,ncol=ny+2)
  x[2:(nx+1),2:(ny+1)]=X     # set the X(data) as a initial 
  x.init=x[2:(nx+1),2:(ny+1)]
  iestop=istop*nx*ny
  for(i1 in 1:iestop) {
    i=sample(1:nx,1)+1
    j=sample(1:ny,1)+1
    # The only other choice for a move is -x, and the proposal distribution
    # is simply the random selection of a new pixel, which is symmetric by
    # default, so the jump probability simplifies to the form below:
    r=exp(-2*beta*x[i,j]*(x[i,j-1]+x[i,j+1]+x[i-1,j]+x[i+1,j]))
    if(runif(1)<r) x[i,j]=-x[i,j]
  }
  gridout <- x[2:(nx+1),2:(ny+1)]
  return(gridout)
}

Energy <- function(Y){            
  m <- dim(Y)[1];n <- dim(Y)[2]
  s1 <- 0;s2 <- 0
  
  for(i in 1:(m-1)){s1 <- s1 +  sum( Y[i,]*Y[(i+1),] ) }
  for(j in 1:(n-1)){s2 <- s2 +  sum( Y[,j]*Y[,(j+1)] ) }
  
  #For 2-dim Ising model 
  #E0 <- sum(sum(Y))     
  #E <- c(E0,s1+s2)              
  E <- s1+s2
  return(E)
}

#-------------------------------------------------------------------------------#

Gen_Ising = function(n){
  S.stat = matrix(NA,nrow=1,ncol=n)
  X.theta = matrix(NA,nrow=1,ncol=n)
  X.dat = matrix(NA,nrow=n,ncol=100)
  
  # theta의 prior를 exp로 주긴했는데 parameter를 모르겠음
  # 적당히 1/0.4406으로 정했음
  
  # Generating x and Sx
  
  for(iter in 1:n){
    theta = rexp(1,1/0.4406)
    x = ising.mh(nx = 10, ny=10, beta = theta, istop = 100) 
    Sx = Energy(x)
    x.vec = as.vector(x)
    
    X.dat[iter,] = x.vec
    S.stat[1,iter] = Sx
    X.theta[1,iter] = theta
  }
  
  return(list("data"=X.dat,"Summary_stat"=S.stat,"Theta"=X.theta))
}

set.seed(4)
train = Gen_Ising(n=1e+06)

write.csv(train$data,"train/train_x.csv",row.names=FALSE)
write.csv(train$Summary_stat,"train/train_y.csv",row.names=FALSE)
write.csv(train$Theta,"train/train_theta.csv",row.names=FALSE)

test = Gen_Ising(n=1e+05)

write.csv(test$data,"test/test_x.csv",row.names=FALSE)
write.csv(test$Summary_stat,"test/test_y.csv",row.names=FALSE)
write.csv(test$Theta,"test/test_theta.csv",row.names=FALSE)



val = Gen_Ising(n=1e+05)

write.csv(val$data,"val/val_x.csv",row.names=FALSE)
write.csv(val$Summary_stat,"val/val_y.csv",row.names=FALSE)
write.csv(val$Theta,"val/val_theta.csv",row.names=FALSE)


# 불확실한 부분 : prior의 parameter, istop?

#-------------------------------------------------------------------------------#