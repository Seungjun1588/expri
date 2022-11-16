setwd("C:/Users/user/desktop/jun/SpatioTemp/hw4")

# Inner sampler


prop_circle = function(){
  # output : vector c()
}

pi_r = function(th1=1.5,th2=10,th3=0.5,R=0,xi,xj,r1,r2){ # R=0 is default in this problem.
  r = sum((xi - xj)^2) # L2 norm.
  
  if( r<= R){
    res = 0
  } else if(r <= r1) {
    res = th1 - ((sqrt(th1)/(th2-R))*(r-th2))^2
  } else { # r > r1
    res = 1 + 1/(th3*(r-r2))^2
  }
  return(res)
}


pi_r(th1=2.5,th2=2,th3=0.5,R=0,xi,xj,r1,r2)



some_fn3 = function(th3=0.5,R=0,r){
  res = 1 + 1/(th3*(r-10))^2
  return(res)
}


some_fn2 = function(th1=1.5,th2=10,R=0,r){
  res = th1 - ((sqrt(th1)/(th2-R))*(r-th2))^2
  return(res)
}


x= seq(0,50,by=0.1)
plot(x,some_fn3(r=x),type="l")
lines(x,some_fn2(r=x),type="l",col="blue")




h_x() = function(xi,xj,lamb,k){
  lamb^(n)*(prod(exp(min(k,sum(log(pi_r(xi,xj,th1=1.5,th2=10,th3=0.5,R=0,r1,r2)))))))
}




d_function(x) = function(){
  return(1/lengnth(x))
}

b_fucntion(A) = function(){
  return(1/A) # what is the diff. btw. A and W?
}

delete_sample = function(){
  
}

Inner_sampler =function(p1=0.3,samples){
  iter= 10
  u = runif(n=1,min=0,max=1)
  
  for(i in 1:iter){
    if(u <p1){
      # add a point
      v1 = runif(n=1,min=0,max=1)
      add_samples =c(samples,prop_circle())
      add_ratio = ((p1-1)*(h_x()*d_function()))/(p1*h_x()*b_function())
      if( v1 < add_ratio){
        samples = add_samples
      }
      
    }
    else{ # u > p1
      # delete a point
      len = length(samples)/2
      v2 = runif(n=1,min=0,max=1)
      del_samples = delete_sample(samples)
      del_ratio = ((p1)*(h_x()*b_function()))/((1-p1)*h_x()*d_function())
      if( v2 < del_ratio){
        samples = del_samples
      }
    }
  }
}

