#include <RcppArmadillo.h>
#include<Rcpp.h>
using namespace Rcpp;
using namespace arma;

double Thetaholding(double& t, double& lambda, double& eta){
  double HRholding; 
  if(fabs(t) < lambda){
    HRholding = 0 ;
  }else{
    HRholding = t/(1+eta);
  }
  return HRholding;
}

double Penalty(double t,double lambda, double eta){
  double p=0;
  if(fabs(t) < lambda/(1+eta))
    p=-0.5*t*t+lambda*fabs(t);
  else
    p=0.5*eta*t*t+0.5*lambda*lambda/(1+eta);
  return p;
}


