#include <RcppArmadillo.h>
#include<Rcpp.h>
#include"touwenjianF.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//'  The package named GIST-fitting(Group Iterative Spectrum Thresholding)
//'  This package is made by Ding Fei, Zhou ya and Chen yingyi 
//' @param y signal response who has been centered
//' @param X1 the design matrix who has been normalized
//' @param lambda is the tunning parameterof the hard-ridge panalty
//' @param Omega is the maximum of interation
//' @param w is the relaxation parameter
//' @return the list of betahat--the estimator of the coefficient,object function--vector of the objection function value with the interation,and interation--the interation when the objection value arrive the optimal value
//' @examples 
//' #Generate the simulation Data
//'set.seed(100)
//'  sigma<-1
//'N<-100
//'y<-c()
//'  for(i in 1:100){
//'    s<-0
//'    nz<-c(0.248,0.25,0.252,0.398,0.4)
//'    Ak<-c(2,4,3,3.5,3)
//'    phik<-c(pi/4,pi/6,pi/3,pi/5,pi/2)
//'    tn<-i
//'    for (fk in nz)
//'    {
//'      en<-rnorm(1,0,sigma)
//'      s=s+Ak*cos(2*pi*fk*tn+phik)+en
//'    }
//'    y[i]<-s
//'  }
//'  x<-c()
//'    for(i in 1:250){
//'      for(j in 1:100){
//'        x[j+(i-1)*100]<-cos(pi*i/250*j)
//'      }
//'    }
//'    for(i in 1:250){
//'      for(j in 1:100){
//'        x[25000+j+(i-1)*100]<-sin(pi*i/250*j)
//'      }
//'    }
//'    X<-matrix(x,100,500)
//'#Normalize the data X and Center the data y,so we get the using Data X1 and y1
//'      X1<-matrix(,100,500)
//'        for(i in 1:500){
//'          X1[,i]=(X[,i]-mean(X[,i]))/sqrt(var(X[,i]))
//'        }
//'        y1 = vector()
//'          for(i in 1:length(y)){
//'            y1[i]= y[i]-1/length(y)*sum(y)
//'          }
//'# Get the Results of the Models Using the Functon in Our Package Icluding the GroupIST and FastGroupIst
//'ResultFastGroupIST =FastGroupIST(X1,y1,Omega = 10000,w = 0.5, epsilon = 0.01,eta =0.1,tau0=30,theta=0.5)
//'ResultFastGroupIST$iteration
//'[1] 16
//'ResultFastGroupIST$`object function`[16]
//'[1] 9.193812
//'sum(ResultFastGroupIST$betahat!=0)
//'[1] 102
//[[Rcpp::export]]
Rcpp::List FastGroupIST(const arma::mat& X1, const arma::vec& y1,double& Omega, double& w, double& epsilon,double& eta,double& tau0 ,double& theta){
  mat X = X1/tau0;
  vec y = y1/tau0;
  int j = 0;
  int n = X.n_cols;
  int N = X.n_rows;
  int i;
  double mj;
  int ite=0;
  //objective function value
  vec fvalue =zeros<vec>(Omega);
  vec beta0 = zeros<vec>(n);
  vec beta1 = zeros<vec>(n);
  vec kesi0 = zeros<vec>(n);
  vec kesi1 = zeros<vec>(n);
  vec l1 = zeros<vec>(n);
  vec middle = zeros<vec>(n);
  double lambda;
  do {
    beta0 = beta1;
    if(j > 0){
      kesi1  = ( 1 - w)*kesi0 + w * (beta0 + X.t()*(y - X * beta0));
    }else{
      kesi1 = beta0 + X.t()*(y - X * beta0);
    }
    mj = N*theta;
    for(i=0 ; i < n/2; i++){
      l1(i) = sqrt(kesi1(i)*kesi1(i)+kesi1(i+n/2)*kesi1(i+n/2));}
    middle =sort(l1,"descend");
    lambda = (middle(mj-1)+middle(mj+1))/2;
    for(i=0 ; i < n/2; i++){  
      if(l1(i) != 0){
        beta1(i) =  kesi1(i) * ThetaholdingF(l1(i),lambda,eta)/l1(i);
        beta1(i+n/2) = kesi1(i+n/2)*ThetaholdingF(l1(i),lambda,eta)/l1(i);
      }else{
        beta1(i)=beta1(i+n/2) = 0;
      }
    }
    //to calculate Penalty term
    double Penaltyterm=0;
    for(i=0 ; i < n/2; i++){
      Penaltyterm=Penaltyterm+PenaltyF(sqrt(beta1(i)*beta1(i)+beta1(i+n/2)*beta1(i+n/2)),lambda,eta);
    }
    fvalue(j)=0.5*norm((tau0*y-tau0*X*beta1),2)+Penaltyterm;
    
    kesi0 = kesi1;
    j = j+1;
    ite=ite+1;
  }
  while ((j < Omega) && (norm(beta1-beta0,2) > epsilon));
  return Rcpp::List::create(Rcpp::Named("betahat")=beta1,
                            Rcpp::Named("object function")=fvalue,
                            Rcpp::Named("iteration")=ite);
  
}
