#include <RcppArmadillo.h>


using namespace Rcpp;
using namespace arma;


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
//init the starting point
mat init(vec size){
  int nrow = size(0), ncol= size(1);
  mat res(nrow,ncol);
  
  for(int i = 0; i< nrow; i++){
    for(int j =0; j< ncol; j++){
      NumericVector num = Rcpp::rnorm( 1,0.0,1.0 );
      IntegerVector res1 = ifelse( num < 0, -1, 1);
      double res2 = as<double>(res1);
      
      res(i,j) = res2 ;
    }
  }
  
  return res;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]  
mat Gen_X(vec size,int n=100){
  //placeholder
  mat X_dat(n,100);


  //Simulation
  double idx = 0;
  for(int i=0; i<n; i++){
    if(i % (n/20) == 0){
      Rcout << (idx/n)*100 << "% done."<< "\n";
    }
    idx += 1;
    
    //init matrix
    mat start = init(size);
    NumericVector beta = rexp(1,0.3);
    double theta = as<double>(beta);
    mat res = MH(start,theta,100);

    X_dat.row(i) = vectorise(res).t();

  }
  return X_dat ;
}

//Rcpp::as<arma::vec>(OUT["vec"])


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]  
mat Gen_theta(vec size,int n=100){
  //placeholder
  vec X_theta(n);

  
  //Simulation
  double idx = 0;
  for(int i=0; i<n; i++){
    if(i % (n/20) == 0){
      Rcout << (idx/n)*100 << "% done."<< "\n";
    }
    idx += 1;
    //init matrix
    mat start = init(size);
    NumericVector beta = rexp(1,0.4406);
    double theta = as<double>(beta);
    mat res = MH(start,theta,100);
    
    X_theta(i) = theta;
  }
  return X_theta ;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]  
mat Gen_S_stat(vec size,int n=100){
  //placeholder
  vec S_stat(n);

  
  //Simulation
  double idx = 0;
  for(int i=0; i<n; i++){
    if(i % (n/20) == 0){
      Rcout << (idx/n)*100 << "% done."<< "\n";
    }
    idx += 1;
    //init matrix
    mat start = init(size);
    NumericVector beta = rexp(1,0.3);
    double theta = as<double>(beta);
    mat res = MH(start,theta,100);
    
    S_stat(i) = Energy(res);
    
  }
  return S_stat ;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List Gen_data(vec size,int n=100){
  //placeholder
  mat X_dat(n,size(0)*size(1));
  vec S_stat(n);
  vec X_theta(n);


  //Simulation
  double idx = 0;
  for(int i=0; i<n; i++){
    if(i % (n/20) == 0){
      Rcout << (idx/n)*100 << "% done."<< "\n";
    }
    idx += 1;
    //init matrix
    mat start = init(size);
    NumericVector beta = rexp(1,0.4);
    double theta = as<double>(beta);
    mat res = MH(start,theta,100);

    X_dat.row(i) = vectorise(res).t();
    X_theta(i) = theta;
    S_stat(i) = Energy(res);

  }
  return List::create(Named("X",X_dat),Named("theta",X_theta),Named("Summary_stat",S_stat)) ;
}





