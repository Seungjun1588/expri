#include <RcppArmadillo.h>
#include <limits>
#include <omp.h>


using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Generating poisson R.V
int rPois(double mu){
  double L = exp(-mu), p = 1;
  int k = 0;
  while( p > L ){
    k = k + 1;
    p = p*randu();
  }
  return(k-1);	
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// calculating poisson density
double dPois(double k, double mu){
  double result, denom = k, K=k;
  
  if( k == 0 ){ result = exp(-mu); }else{
    while( k > 1  ){
      k = k-1;
      denom = denom*(k);	
    }
    result = pow(mu,K)*exp(-mu)/denom;
  }
  return(result);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Sign fucntion
int Sign(double x){
  int result;
  if( x < 0 ){ result = -1; }else{
    if( x > 0 ){ result = 1; }else{
      result = 0;
    }}
  
  return(result);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Summary statistics for Ising model
int Energy(mat X){
  int nrow = X.n_rows, ncol = X.n_cols, s1 = 0, s2 = 0;
  int result;
  
  for(int i = 0; i< nrow-1; i++){
    s1 = s1 + accu( X.row(i)%X.row(i+1) );
  }
  for(int j = 0; j< ncol-1; j++){
    s2 = s2 + accu( X.col(j)%X.col(j+1) );
  }
  
  result = s1 + s2;
  
  return(result);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Perfect sampling
mat ProppWilson(int nrow, int ncol, double b){
  mat Xmin = -1*ones(nrow,ncol), Xmax = ones(nrow,ncol);
  mat work1 = Xmin, work2 = Xmax;
  work1.insert_cols(0,zeros(nrow));
  work1.insert_cols(ncol+1,zeros(nrow));
  work1.insert_rows(0,trans(zeros(ncol+2)));
  work1.insert_rows(nrow+1,trans(zeros(ncol+2)));
  work2.insert_cols(0,zeros(nrow));
  work2.insert_cols(ncol+1,zeros(nrow));
  work2.insert_rows(0,trans(zeros(ncol+2)));
  work2.insert_rows(nrow+1,trans(zeros(ncol+2)));	
  int Time = 1; 
  
  // until Xmin == Xmax  
  while(  accu(abs(Xmax-Xmin)) > 0  ){    	
    Time = 2*Time;
    for(int k = 0; k< Time; k++){
      // just one component random scan update
      int i = ceil(nrow*randu());   
      int j = ceil(ncol*randu());      
      double u  = randu();		
      // update for Xmin
      double r1 = exp( 2*b*( work1(i,j-1)+work1(i,j+1)+work1(i-1,j)+work1(i+1,j) ) );
      double p1 = r1/(1+r1);                
      if( u < p1  ){
        Xmin( (i-1),(j-1) ) = 1;
        work1(i,j) = 1;
      }else{
        Xmin( (i-1),(j-1) ) = -1;	
        work1(i,j) = -1;
      }
      // update for Xmax
      double r2 = exp( 2*b*( work2(i,j-1)+work2(i,j+1)+work2(i-1,j)+work2(i+1,j) ) );
      double p2 = r2/(1+r2);  
      if( u < p2  ){
        Xmax( (i-1),(j-1) ) = 1;
        work2(i,j) = 1;
      }else{
        Xmax( (i-1),(j-1) ) = -1;	
        work2(i,j) = -1;
      }
    } 
    
  }
  
  return(Xmin);
} 


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// random scan gibbs update
mat Gibb(mat initial, double b, double cycle){
  int nrow = initial.n_rows, ncol = initial.n_cols;
  mat work = initial;			
  work.insert_cols(0,zeros(nrow));
  work.insert_cols(ncol+1,zeros(nrow));
  work.insert_rows(0,trans(zeros(ncol+2)));
  work.insert_rows(nrow+1,trans(zeros(ncol+2)));
  int iestop = cycle*nrow*ncol;
  
  for(int k = 0; k< iestop; k++){
    
    int i = ceil(nrow*randu());   
    int j = ceil(ncol*randu());  
    
    double r = exp( 2*b*( work(i,j-1)+work(i,j+1)+work(i-1,j)+work(i+1,j) ) );
    double p = r/(1+r);                  
    if( randu() < p  ){
      initial( (i-1),(j-1) ) = 1;
      work(i,j) = 1;
    }else{
      initial( (i-1),(j-1) ) = -1;	
      work(i,j) = -1;
    }
  }
  
  return(initial);	
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// MH random scan update 
mat MH(mat initial, double b, double cycle){
  int nrow = initial.n_rows, ncol = initial.n_cols;
  mat work = initial;
  work.insert_cols(0,zeros(nrow));
  work.insert_cols(ncol+1,zeros(nrow));
  work.insert_rows(0,trans(zeros(ncol+2)));
  work.insert_rows(nrow+1,trans(zeros(ncol+2)));
  int iestop = cycle*nrow*ncol;
  
  for(int k = 0; k< iestop; k++){
    int i = ceil(nrow*randu());   
    int j = ceil(ncol*randu());  
    double p = exp( -2*b*work(i,j)*(work(i,j-1)+work(i,j+1)+work(i-1,j)+work(i+1,j)) );
    
    if( randu() < p  ){
      initial( (i-1),(j-1) ) = -initial( (i-1),(j-1) );
      work(i,j) = -work(i,j) ;
    }else{
      initial( (i-1),(j-1) ) = initial( (i-1),(j-1) );	
      work(i,j) = work(i,j);
    }	
  }  
  
  return(initial);	
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Pseudo likelihood mcmc for Ising model example
vec IsingPseudoMCMC(int outer, double initial, double sigma, mat X){
  vec parameter(outer);
  parameter(0) = initial;
  double logprob,u,bprop;
  double negativeInf = -std::numeric_limits<float>::infinity();;	 
  int nrow = X.n_rows, ncol = X.n_cols;
  mat work = X;
  work.insert_cols(0,zeros(nrow));
  work.insert_cols(ncol+1,zeros(nrow));
  work.insert_rows(0,trans(zeros(ncol+2)));
  work.insert_rows(nrow+1,trans(zeros(ncol+2)));
  
  double loglike=0, loglikeproposed=0;
  for(int i = 1; i< nrow+1; i++){
    for(int j = 1; j< ncol+1; j++){
      double p = exp( 2*parameter(0)*(work(i,j-1)+work(i,j+1)+work(i-1,j)+work(i+1,j)) );
      p = p/(1+p);
      loglike = loglike + ((work(i,j)+1)/2)*log(p) + (1-(work(i,j)+1)/2)*log(1-p);
    }
  }
  
  // Start of MCMC Chain 
  for(int k = 0; k< outer-1; k++){
    // propose parameters 
    bprop = parameter(k) + sigma*randn();
    
    if( bprop > 1 || bprop < 0 ){
      logprob = negativeInf;	
    }else{
      loglikeproposed = 0; 
      for(int i = 1; i< nrow+1; i++){
        for(int j = 1; j< ncol+1; j++){
          double p = exp( 2*bprop*(work(i,j-1)+work(i,j+1)+work(i-1,j)+work(i+1,j)) );
          p = p/(1+p);
          loglikeproposed = loglikeproposed + ((work(i,j)+1)/2)*log(p) + (1-(work(i,j)+1)/2)*log(1-p);
        }
      } 
      logprob = loglikeproposed - loglike;  
    }
    u = log( randu() );
    if( u< logprob ){
      parameter(k+1) = bprop;
      loglike = loglikeproposed;
    }else{
      parameter(k+1) = parameter(k);
      loglike = loglike;
    } 
    
  }  
  
  return(parameter);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Moller's auxiliary variable algorithm for Ising model
vec IsingMoller(int outer, double initial, double sigma, mat X){
  vec parameter(outer);
  parameter(0) = initial;
  double logprob,u,bprop;
  double negativeInf = -std::numeric_limits<float>::infinity();;	
  int nrow = X.n_rows, ncol = X.n_cols;
  
  // auxiliary variable 
  mat Y = ProppWilson(nrow,ncol,initial); 
  int XStat = Energy(X), YStat = Energy(Y);
  
  // Start of MCMC Chain 
  for(int k = 0; k< outer-1; k++){
    // propose parameters 
    bprop = parameter(k) + sigma*randn();
    // propose auxiliary variables 
    mat Yprop = ProppWilson(nrow,ncol,bprop);
    int YpropStat = Energy(Yprop);
    
    if( bprop > 1 || bprop < 0 ){
      logprob = negativeInf;	
    }else{
      
      logprob = 
        ( initial*YpropStat + bprop*XStat + parameter(k)*YStat ) -
        ( initial*YStat + parameter(k)*XStat + bprop*YpropStat );
    }
    
    u = log( randu() );
    if( u< logprob ){
      parameter(k+1) = bprop;
      YStat = YpropStat;
    }else{
      parameter(k+1) = parameter(k);
      YStat = YStat;
    } 
    
  }  
  
  return(parameter);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Exchange algorithm for Ising model
vec IsingExchange(int outer, double initial,double sigma, mat X){
  vec parameter(outer);
  parameter(0) = initial;
  double logprob,u,bprop;
  double negativeInf = -std::numeric_limits<float>::infinity();;	
  int nrow = X.n_rows, ncol = X.n_cols;
  int XStat = Energy(X), YStat;
  mat Y(nrow,ncol); // auxiliary variable
  
  // Start of MCMC Chain 
  for(int k = 0; k< outer-1; k++){
    // propose parameters 
    bprop = parameter(k) + sigma*randn();
    // propose auxiliary variables 
    Y = ProppWilson(nrow,ncol,bprop);
    YStat = Energy(Y);
    
    if( bprop > 1 || bprop < 0 ){
      logprob = negativeInf;	
    }else{
      logprob = ( bprop-parameter(k) )*XStat +   ( parameter(k)-bprop )*YStat;
    }
    u = log( randu() );
    if( u< logprob ){
      parameter(k+1) = bprop;
    }else{
      parameter(k+1) = parameter(k);
    } 
    
  }  
  
  return(parameter);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Double Metropolis Hastings for Ising model
vec IsingDMH(int outer, int inner, double initial, double sigma, mat X){
  vec parameter(outer);
  parameter(0) = initial;
  double logprob,u,bprop;
  double negativeInf = -std::numeric_limits<float>::infinity();;	
  int nrow = X.n_rows, ncol = X.n_cols;
  int XStat = Energy(X), YStat;
  mat Y(nrow,ncol); // auxiliary variable
  
  // Start of MCMC Chain 
  for(int k = 0; k< outer-1; k++){
    // propose parameters 
    bprop = parameter(k) + sigma*randn();
    // propose auxiliary variables 
    Y = Gibb(X,bprop,inner);
    YStat = Energy(Y);
    
    if( bprop > 1 || bprop < 0 ){
      logprob = negativeInf;	
    }else{
      logprob = ( bprop-parameter(k) )*XStat + ( parameter(k)-bprop )*YStat;
    }
    u = log( randu() );
    if( u< logprob ){
      parameter(k+1) = bprop;
    }else{
      parameter(k+1) = parameter(k);
    } 
    
  }  
  
  return(parameter);
}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Atchade's adpative algorithm for Ising model
vec IsingAtchade(int outer, int inner, double initial, double sigma, vec th, mat X){	    
  int d = th.n_rows, p = 1, nrow = X.n_rows, ncol = X.n_cols, Vmin;
  vec Vis= zeros(d), gamVec(4), cVec;    // gamVec length and sequence both can be changed 
  for(int i = 0; i< 4; i++){ gamVec[i] = pow( 0.1, i ); }
  
  vec parameter(outer); 
  mat Data = X;
  int Stat0 = Energy(X), Stat;
  parameter(0) = initial;
  mat Esum = zeros(outer,d); // summary statistics will be stored
  
  // approximate cVec ( logZ(theta_{i}) ) until gam is small
  for(int igam = 0; igam< 4; igam++){  
    
    if( igam == 3 ){ Vmin=1000; }else{ Vmin=0; } //Vmin=1000 can be changed
    double gam = gamVec(igam);
    X = Data, Stat = Stat0; 
    vec VisTemp = zeros(d);     
    
    // pos denotes the variable I in Algorithm 2.1
    int pos = d-1;                          
    
    // cVec is the variable c in Algorithm 2.1
    if( igam == 0 ){ cVec = zeros(d); }
    
    // Stopping criteria
    while( VisTemp.min() <= Vmin || abs(  VisTemp-mean(VisTemp)*ones(d)  ).max() >= 0.2*mean(VisTemp) ){
      // Update X_{n}
      X = MH(X, th(pos), inner);
      Stat = Energy(X);
      
      // Update I; meaning pos
      double mx = (th*Stat - cVec).max();
      vec A = exp((th*Stat - cVec - mx)) / sum( exp((th*Stat - cVec - mx))	);
      double u = randu(), l = -1, om = 0;
      while( om<u ){
        l = l+1;
        om = om +A[l];
      }
      pos = l;
      
      // Update c
      cVec = cVec + log(1+gam)*A;
      
      // We need the next two updates to control the sequence gamma
      VisTemp[pos] = VisTemp[pos] + 1; 		
      //if( igam ==3 ){		
      //Vis[pos] = Vis[pos] + 1;
      //Esum(Vis[pos]-1,pos) = Stat; 			
      //}
      
    }
  }  
  
  // We are now ready the start the MCMC
  // From here, gamma is deterministic small value (0.001) 
  // Initilization of ${theta_n},X,I,c$    
  double theta, thetaprop;
  theta =  parameter(0);
  X = Data, Stat = Stat0;
  int pos = d-1;    
  double gam = gamVec(3),prob; 
  
  // Proposal Parameters
  double b = -1.5, Acc = 0, tau = 0.3;
  
  // c = cVec from above iteration
  // Bandwidth for smoothing the estimate
  int hpar = 10;
  for(int k = 0; k< outer-1; k++){
    
    // Update X_{n}
    X = MH(X, th(pos), inner);
    Stat = Energy(X);
    
    // Update I; meaning pos
    double mx = (Stat*th - cVec).max();
    vec A = exp((Stat*th - cVec - mx)) / sum( exp((Stat*th - cVec - mx))	);
    double u = randu(), l = -1, om = 0;
    while( om<u ){
      l = l+1;
      om = om +A[l];
    }
    pos = l;
    
    // Update c
    cVec = cVec + log(1+gam)*A;
    Vis[pos] = Vis[pos] + 1;
    Esum(Vis[pos]-1,pos) = Stat; 
    
    // Update theta if there are enough observations
    // Evaluate \pi(theta)
    // Distance to each theta^i
    uvec ind = sort_index(  abs( theta*ones(d) - th )    );   	
    vec cw = (1/hpar)*ones(hpar);
    
    // There might be some trouble here if none of the closest 'hpar'
    // particles around thetaprop has received data
    vec expEtheta = zeros(hpar);
    for(int i = 0; i< hpar; i++){	
      if( Vis[ind[i]] == 0 ){ expEtheta[i] = 0; }else{     
        vec tmp = ( theta - th[ind[i]] )*Esum.col( ind[i] ) ;	
        tmp = tmp(  span( 0, Vis[ind[i]]-1 )  );
        double tmpMax = tmp.max();	
        expEtheta[i] = tmpMax + log( sum( exp( tmp - tmpMax ) ) ) - log( Vis[ind[i]] ); // eq(8)'s [   ] part       	
      }
    }
    vec dummy2(hpar);
    for(int i = 0; i< hpar; i++){ dummy2[i] = cVec[ ind( span(0,hpar-1) )[i] ] + expEtheta[i] + cw[i];	}
    double ftheta = dummy2.max();                        
    // Eth = log(exp(E(x,theta))) - log(Z(theta)) 
    // Stat0*theta: log(exp(E(x,theta)))= E(x,theta)
    // log(sum( exp(cVec[ind[1:hpar]]+expEtheta+cw) )): log(Z(theta))  
    double Eth = Stat0*theta - ftheta - log(sum( exp(dummy2 - ftheta) ));
    
    // Propose a new theta
    thetaprop = theta + sigma*randn();
    if( thetaprop > 1 || thetaprop < 0 ){ prob = 0; }else{
      
      // Evaluate posterior at the new theta        
      ind = sort_index(  abs( thetaprop*ones(d) - th )    );   	
      cw = (1/hpar)*ones(hpar);
      
      // There might be some trouble here if none of the closest 'hpar'
      // particles around thetaprop has received data
      vec expEthetaprop = zeros(hpar);        
      for(int i = 0; i< hpar; i++){	
        if( Vis[ind[i]] == 0 ){ expEthetaprop[i] = 0; }else{
          vec tmp = ( thetaprop - th[ind[i]] )*Esum.col( ind[i] ) ;	
          tmp = tmp(  span( 0, Vis[ind[i]]-1 )  );
          double tmpMax = tmp.max();	
          expEthetaprop[i] = tmpMax + log( sum( exp( tmp - tmpMax ) ) ) - log( Vis[ind[i]] ); // eq(8)'s [   ] part       	
        }
      }
      for(int i = 0; i< hpar; i++){ dummy2[i] = cVec[ ind( span(0,hpar-1) )[i] ] + expEthetaprop[i] + cw[i];	}
      double fthetaprop = dummy2.max();   
      double Ethprop = Stat0*thetaprop - fthetaprop - log(sum( exp(dummy2 - fthetaprop) ));						
      
      // Acceptance prob.
      prob = exp(Ethprop-Eth);
    }
    // Accept - Reject
    u = randu();
    if(u <= prob){ theta = thetaprop; }
    
    // Adaptive scaling of the proposal
    vec MIN(2);
    MIN[0] = 1, MIN[1] = prob;
    Acc = Acc + (1/(k+1))*( MIN.min() - Acc ); 
    b = b +(1/(k+1))*( MIN.min() - tau );
    
    // Save output
    parameter(k+1) = theta;
  }
  
  
  
  return(parameter);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Anealed importance sampling logZ estimate for nrow by ncol Ising model
double logZAIS(int noImp, int noTemp, double b, int nrow, int ncol){
  // noImp - number of importance samples
  // a and b - Ising params
  // importance samples done in parallel
  
  // temperature schedule
  vec beta = linspace<vec>(0, 1, noTemp+1);
  vec logweight(noImp);
  
  int M;
#pragma omp parallel shared(logweight) private(M)
{	
#pragma omp for schedule(static)  
  for(M = 0; M< noImp; M++){
    
    // sample random grid 
    mat y = sign(randn(nrow,ncol));
    mat work = y;			
    work.insert_cols(0,zeros(nrow));
    work.insert_cols(ncol+1,zeros(nrow));
    work.insert_rows(0,trans(zeros(ncol+2)));
    work.insert_rows(nrow+1,trans(zeros(ncol+2)));
    
    int Stat = Energy(y);
    double logw = (beta[1]-beta[0])*b*Stat - (beta[0]-beta[1])*(nrow*ncol)*log(2);
    
    // sample from succesive tempered distributions and calc weight
    for(int n = 1; n< noTemp; n++){
      // random scan Gibbs sampler
      // choose one spin in each grid to update
      int i = ceil(nrow*randu());   
      int j = ceil(ncol*randu());   	   
      double r = exp( 2*beta[n]*b*( work(i,j-1)+work(i,j+1)+work(i-1,j)+work(i+1,j) ) );
      double p = r/(1+r); 
      
      if( randu() < p  ){
        y( (i-1),(j-1) ) = 1;
        work(i,j) = 1;
      }else{
        y( (i-1),(j-1) ) = -1;	
        work(i,j) = -1;
      }
      // calc AIS weight in log space
      Stat = Energy(y);
      logw = logw + (beta[n+1]-beta[n])*b*Stat - (beta[n]-beta[n+1])*(nrow*ncol)*log(2);
    }
    logweight[M] = logw;     
  }
}

double mx = logweight.max();
double logZ = log(mean(exp(logweight - mx))) + mx; 	
return(logZ);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Playing Russian Roulette
List Roulette(int nrow, int ncol, double logZtilde, int NimpSamps, double bprop, double r, double cx, int cmax, int ru, double mu, int NTemp){
  double series, term, weight, maxk, Corr;
  int count, k;
  vec logZx, kappa, q;
  
  if( ru == 1 ){
    // use roulette truncation
    while( ru == 1 ){
      count = 1, series = 1, term = 1, weight = 1;
      logZx = zeros(cmax), kappa = zeros(cmax), q = zeros(cmax);
      while( count < cmax ){
        
        // calculate estimate of Z and kappa
        logZx[count-1] = logZAIS(NimpSamps, NTemp, bprop, nrow,ncol);
        kappa[count-1] = 1  -  cx * exp(logZx[count-1]-logZtilde);
        
        // if kappa is outside convergence region make c and r smaller
        if( fabs(kappa[count-1]) > 1 ){ break; }			  	
        term = term*kappa[count-1];
        
        if( fabs(term) < r ){
          // play Roulette
          q[count-1] = fabs(term)/r;
          weight = weight*q[count-1];    
          if( q[count-1] > randu() ){
            // survivie
            count = count + 1;
            series = series + (term/weight);	
          }else{
            // terminate
            maxk = count - 1; 
            ru = 0;
            break;}}else{
              // add to go to the next term
              series = series + term;
              count = count + 1; 
            }
            
      } // end of while (count < cmax)
      
      // leave while loop if cmax exceeded or roulette terminated
      if( count >= cmax || ru == 0 ){ break; }else{
        // otherwise make c smaller to ensure convergence    
        cx = cx*0.9;
        r = 1-cx;
        // In R, remove logZx, kappa, q at here
      }  
    } // end of while(ru == 1)   
    
    // if roulette terminated in time
    if( count < cmax ){ Corr = series; }else{	
      maxk = cmax;
      // use poisson truncation, Corr is calculated through inverse of poisson density;
      k = rPois(mu);
      if( k == 0 ){ Corr = 1/dPois(k,mu); }else{
        while( ru == 0 ){
          kappa = zeros(k), logZx = zeros(k);
          for(int ii = 0; ii < k; ii++){
            logZx[ii] = logZAIS(NimpSamps, NTemp, bprop, nrow,ncol);
            kappa[ii] = 1  -  cx * exp(logZx[ii]-logZtilde);     	 	   	
            if( ( ii == k) && ( fabs(kappa[ii]) < 1 ) ){ ru = 1; } 
            if( fabs(kappa[ii]) > 1 ){
              cx = cx*0.9;
              // In R, remove logZx, kappa at here
              break;
            }
          }
        } // end of while( ru == 0 )
        
        Corr = prod(kappa)/dPois(k,mu);  
      }  
    }
  }else{
    // use poisson truncation
    k = rPois(mu);
    maxk = k;
    if( k == 0 ){ Corr = 1/dPois(k,mu); }else{      	
      while( ru == 0 ){
        kappa = zeros(k), logZx = zeros(k); 
        for(int ii = 0; ii < k; ii++){       
          logZx[ii] = logZAIS(NimpSamps, NTemp, bprop, nrow,ncol);
          kappa[ii] = 1  -  cx * exp(logZx[ii]-logZtilde);     	 
          if( ( ii == k) && ( fabs(kappa[ii]) < 1 ) ){ ru = 1; } 
          if( fabs(kappa[ii]) > 1 ){
            cx = cx*0.9;
            // In R, remove logZx, kappa at here
            break;
          }
        }
      } // end of while( ru == 0 )
      
      Corr = prod(kappa)/dPois(k,mu);  
    } 
  }
  
  return List::create(Named("Corr") = Corr, Named("c") = cx, Named("maxk") = maxk);	
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Russian Roulette for Ising model
mat IsingRoulette(int outer, int NimpSamps, int Ntemp, double initial, double sigma, double r, double cx, int cmax, int ru, double mu, mat X){
  int nrow = X.n_rows, ncol = X.n_cols, Stat = Energy(X), sCorrprop;
  double lhX, lhXp, logZtilde, logprob,u,bprop, Corr, c, maxk;
  double negativeInf = -std::numeric_limits<float>::infinity();;;
  
  
  // posterior sample, sign will be stored. 
  mat parameter(2,outer); 
  // initial values
  parameter(0,0) = initial, parameter(1,0) = 1;
  
  // current log likelihood
  logZtilde = logZAIS(NimpSamps, Ntemp, parameter(0,0), nrow, ncol);
  lhX = parameter(0,0)*Stat + log(cx) - logZtilde;
  
  // MCMC chain 
  for(int k = 0; k< outer-1; k++){
    
    // propose parameters 
    bprop = parameter(0,k) + sigma*randn();  
    if( bprop > 1 || bprop < 0 ){ logprob = negativeInf; }else{
      
      // calculate log Z_tilde using extra importance samples
      logZtilde = logZAIS(2*NimpSamps, Ntemp, bprop, nrow, ncol);
      
      // calculate correction  
      List output = Roulette(nrow, ncol, logZtilde, NimpSamps, bprop, r, cx, cmax, ru, mu, Ntemp);
      vec output1 = output[0], output2 = output[1], output3 = output[2];
      Corr = output1[0], c = output2[0], maxk = output3[0];
      sCorrprop = Sign(Corr);
      
      // calculate estimated likelihood in log-space
      lhXp = bprop*Stat + log(c) + log(fabs(Corr)) - logZtilde;
      // calculate MH ratio
      logprob = lhXp - lhX;
    }
    u = log( randu() );
    if( u < logprob ){
      parameter(0,k+1) = bprop;
      parameter(1,k+1) = sCorrprop;
      lhX = lhXp;
    }else{
      parameter(0,k+1) = parameter(0,k);
      parameter(1,k+1) = parameter(1,k);
      lhX = lhX;
    }
  }
  
  
  return(parameter);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Noisy Exchange algorithm for Ising model
vec IsingNExchange(int outer, int NimpSamps, double initial, double sigma, mat X){
  vec parameter(outer), logweight(NimpSamps);
  parameter(0) = initial;
  double logprob,u,bprop, mx, logRatio;
  double negativeInf = -std::numeric_limits<float>::infinity();;	
  int nrow = X.n_rows, ncol = X.n_cols;
  int XStat = Energy(X);
  
  
  // Start of MCMC chain 
  for(int k = 0; k< outer-1; k++){
    // propose parameters 
    bprop = parameter(k) + sigma*randn();
    if( bprop > 1 || bprop < 0 ){ logprob = negativeInf; }else{
      
      int M;
#pragma omp parallel shared(logweight) private(M)
{	
#pragma omp for schedule(static)  
  for(M = 0; M< NimpSamps; M++){
    mat Y = ProppWilson(nrow,ncol,bprop); // auxiliary variable
    logweight[M] = ( parameter(k)-bprop )*Energy(Y);
  }
} // End of parallelization       

mx = logweight.max();
logRatio = log(mean(exp(logweight - mx))) + mx; 	
logprob = ( bprop-parameter(k) )*XStat + logRatio;
    }
    u = log( randu() );
    if( u< logprob ){
      parameter(k+1) = bprop;
    }else{
      parameter(k+1) = parameter(k);
    }
    
  }
  
  return(parameter);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Noisy DMH algorithm for Ising model
vec IsingNDMH(int outer, int inner, int NimpSamps, double initial, double sigma, mat X){
  vec parameter(outer), logweight(NimpSamps);
  parameter(0) = initial;
  double logprob,u,bprop, mx, logRatio;
  double negativeInf = -std::numeric_limits<float>::infinity();;	
  int nrow = X.n_rows, ncol = X.n_cols;
  int XStat = Energy(X);
  
  
  // Start of MCMC chain 
  for(int k = 0; k< outer-1; k++){
    // propose parameters 
    bprop = parameter(k) + sigma*randn();
    if( bprop > 1 || bprop < 0 ){ logprob = negativeInf; }else{
      
      int M;
#pragma omp parallel shared(logweight) private(M)
{	
#pragma omp for schedule(static)  
  for(M = 0; M< NimpSamps; M++){
    mat Y = Gibb(X,bprop,inner);  // auxiliary variable
    logweight[M] = ( parameter(k)-bprop )*Energy(Y);
  }
} // End of parallelization       

mx = logweight.max();
logRatio = log(mean(exp(logweight - mx))) + mx; 	
logprob = ( bprop-parameter(k) )*XStat + logRatio;
    }
    u = log( randu() );
    if( u< logprob ){
      parameter(k+1) = bprop;
    }else{
      parameter(k+1) = parameter(k);
    }
    
  }
  
  return(parameter);
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Count visitation for Atchade's adpative algorithm
vec StopAtchade(int inner, vec th, mat X){	    
  int d = th.n_rows, p = 1, nrow = X.n_rows, ncol = X.n_cols, Vmin;
  vec Vis= zeros(d), VisTemp = zeros(d), gamVec(4), cVec;    // gamVec length and sequence both can be changed 
  for(int i = 0; i< 4; i++){ gamVec[i] = pow( 0.1, i ); }
  
  mat Data = X;
  int Stat0 = Energy(X), Stat;
  
  // approximate cVec ( logZ(theta_{i}) ) until gam is small
  for(int igam = 0; igam< 4; igam++){  
    
    if( igam == 3 ){ Vmin=1000; }else{ Vmin=0; } //Vmin=1000 can be changed
    double gam = gamVec(igam);
    X = Data, Stat = Stat0; 
    VisTemp = zeros(d);     
    
    // pos denotes the variable I in Algorithm 2.1
    int pos = d-1;                          
    
    // cVec is the variable c in Algorithm 2.1
    if( igam == 0 ){ cVec = zeros(d); }
    
    // Stopping criteria
    while( VisTemp.min() <= Vmin || abs(  VisTemp-mean(VisTemp)*ones(d)  ).max() >= 0.2*mean(VisTemp) ){
      // Update X_{n}
      X = MH(X, th(pos), inner);
      Stat = Energy(X);
      
      // Update I; meaning pos
      double mx = (th*Stat - cVec).max();
      vec A = exp((th*Stat - cVec - mx)) / sum( exp((th*Stat - cVec - mx))	);
      double u = randu(), l = -1, om = 0;
      while( om<u ){
        l = l+1;
        om = om +A[l];
      }
      pos = l;
      
      // Update c
      cVec = cVec + log(1+gam)*A;	
      
      // We need the next two updates to control the sequence gamma
      VisTemp[pos] = VisTemp[pos] + 1; 		
      //if( igam ==3 ){		
      //Vis[pos] = Vis[pos] + 1;
      //Esum(Vis[pos]-1,pos) = Stat; 			
      //}
      
    }
    Vis = Vis + VisTemp;
  }
  return(Vis);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
vec IsingFDMH(int outer, int inner, double initial, double sigma, mat X){
  vec parameter(outer);
  parameter(0) = initial;
  double logprob,u,bprop;
  double negativeInf = -std::numeric_limits<float>::infinity();; 
  int nrow = X.n_rows, ncol = X.n_cols;
  int XStat = Energy(X), YStat;
  mat Y(nrow,ncol); // auxiliary variable
  
  // Start of MCMC Chain 
  for(int k = 0; k< outer-1; k++){
    // propose parameters 
    bprop = parameter(k) + sigma*randn();
    // propose auxiliary variables 
    Y = Gibb(X,bprop,inner);
    YStat = Energy(Y);
    
    if( bprop > 1 || bprop < 0 ){
      logprob = negativeInf; 
    }else{
      logprob = ( bprop-parameter(k) )*XStat + ( parameter(k)-bprop )*YStat;
      logprob = 0.5*logprob;
    }
    u = log( randu() );
    if( u< logprob ){
      parameter(k+1) = bprop;
    }else{
      parameter(k+1) = parameter(k);
    }    
  }  
  
  return(parameter);
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// AEX for Ising model
List IsingAEX(int Niter, int Numaux, double cycle, double t0, int neighbor, vec th, mat thdist, double initial, double sigma, mat X){ 
  
  int auxNiter = Numaux*20 + Numaux;
  double d = th.n_rows;                                        // number of particles
  vec p = (1/d)*ones(d), lw = zeros(d), Vis = zeros(d);        // target prob, logartihm of abundance factor
  mat  Sdummy = zeros(auxNiter,3);                                  // data base suff, parameter, abundance factor will be stored
  
  // initialize
  int indprop, ind = 0;
  Sdummy(0,2) = lw[ind], Sdummy(0,1) = th[ind];
  mat auxvar = Gibb(X, Sdummy(0,1), cycle);
  Sdummy(0,0) = Energy( auxvar );
  
  // preliminary run of the auxiliary chain 
  for(int k = 0; k< auxNiter-1; k++){
    // decide update theta or aux var 
    if( randu() < 0.75 ){  // update theta
      uvec dummy = sort_index(  thdist.row( ind )   );  // using distance matrix, ordering from current theta 
      uvec samplist = dummy( span( 1, neighbor )  );    // neighborhood calcul
      int ii  = floor( (neighbor)*randu() ) ;         // sample uniformly from neighborhood  
      indprop = samplist[ii];
      
      double thprop = th[indprop];
      Sdummy(k+1,0) = Sdummy(k,0);
      double logprob = lw[ind] - lw[indprop]  + Sdummy(k,0)*( thprop-Sdummy(k,1) ); 
      
      double u = log( randu() );
      if( u< logprob ){
        Sdummy(k+1,1) = thprop;
        ind = indprop;
      }else{
        Sdummy(k+1,1) = Sdummy(k,1);
        ind = ind; 
      }  
      
    }else{                 // update aux var, always accept by DBE
      auxvar = Gibb(auxvar, Sdummy(k,1), cycle);
      ind = ind, Sdummy(k+1,1) = Sdummy(k,1); Sdummy(k+1,0) = Energy(auxvar);
    }
    
    // update abundance factor
    vec e = zeros(d);
    e[ind] = 1;
    double val;
    double kdummy = k;
    if( t0 > kdummy ){ val = t0; }else{ val = kdummy; }
    lw = lw + (t0/val)*(e-p);
    Sdummy(k+1,2) = lw[ind];
    
    // record visitation frequencies
    Vis = Vis + e;
  }  
  
  // burn in Numaux number and equally 20 sample spaced
  mat S = zeros(Numaux+Niter,3); 
  for(int k = 0; k< Numaux; k++){
    S(k,0) = Sdummy(Numaux-1 + (k+1)*20, 0);   
    S(k,1) = Sdummy(Numaux-1 + (k+1)*20, 1);
    S(k,2) = Sdummy(Numaux-1 + (k+1)*20, 2);
  }
  
  
  // initialize final chain
  vec parameter(Niter);
  parameter(0) = initial;
  double negativeInf = -std::numeric_limits<float>::infinity();; 
  int XStat = Energy(X);
  double logprob;
  
  // run auxiliary chain and target chain simultaneously
  for(int k = 0; k< Niter-1; k++){
    
    // Auxiliary chain 
    // decide update theta or aux var 
    if( randu() < 0.75 ){  // update theta
      uvec dummy = sort_index(  thdist.row( ind )   );  // using distance matrix, ordering from current theta 
      uvec samplist = dummy( span( 1, neighbor )  );    // neighborhood calcul
      int ii  = floor( (neighbor)*randu() ) ;         // sample uniformly from neighborhood  
      indprop = samplist[ii];
      
      double thprop = th[indprop];
      S(Numaux-1 + k + 1, 0) = S(Numaux-1 + k,0);
      double logprob1 = lw[ind] - lw[indprop]  + S(Numaux-1 + k,0)*( thprop-S(Numaux-1 + k,1) ); 
      
      double u = log( randu() );
      if( u< logprob1 ){
        S(Numaux-1 + k + 1, 1) = thprop;
        ind = indprop;
      }else{
        S(Numaux-1 + k + 1, 1) = S(Numaux-1 + k, 1);
        ind = ind; 
      } 
    }else{                 // update aux var, always accept by DBE
      auxvar = Gibb(auxvar, S(Numaux-1 + k, 1), cycle);
      ind = ind, S(Numaux-1 + k + 1, 1) = S(Numaux-1 + k, 1); S(Numaux-1 + k + 1, 0) = Energy(auxvar);
    }
    
    // update abundance factor
    vec e = zeros(d);
    e[ind] = 1;
    double val;
    double kdummy = auxNiter-1 + k;
    if( t0 > kdummy ){ val = t0; }else{ val = kdummy; }
    lw = lw + (t0/val)*(e-p);
    S(Numaux-1 + k + 1, 2) = lw[ind];
    
    
    // Target chain
    double bprop = parameter(k) + sigma*randn();   
    if( bprop > 1 || bprop < 0 ){
      logprob = negativeInf; 
    }else{
      
      vec Prob = zeros(Numaux + k + 1);  // because the number of row of S is Numaux + k + 1
      for(int m = 0; m< Numaux + k + 1; m++){
        Prob[m] =  S(m,2) + S(m,0)*(bprop-S(m,1)); 
      }
      
      double mx = (  Prob  ).max();
      Prob = exp(Prob - mx);                    
      Prob = Prob/sum(Prob);     
      double uu = randu(), l = -1, om = 0;
      while( om<uu ){
        l = l+1;
        om = om + Prob[l];
      } 
      int Statprop = S(l,0);
      logprob = ( bprop-parameter(k) )*XStat + ( parameter(k)-bprop )*Statprop;
    }
    
    if( log( randu() )< logprob ){
      parameter(k+1) = bprop;
    }else{
      parameter(k+1) = parameter(k);
    }     
  }  
  
  
  return List::create(Named("Vis") = Vis, Named("S") = S, Named("abundance") = lw, Named("par") = parameter); 
}

