
#include <RcppEigen.h>
#include <Rcpp.h>
#include <iostream>
#include <algorithm>
// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
Eigen::MatrixXd sigmak(Eigen::MatrixXd Y, Eigen::MatrixXd Y_mean, int k ,int n){
  double nn = n;
  Y = Y-Y_mean.replicate(1,n);
  //Rcout<<Y_mean.replicate(1,n);
  Eigen::MatrixXd Cov_lagk = Y.rightCols(n-k)*Y.leftCols(n-k).transpose()/nn;
  return Cov_lagk;
}
// [[Rcpp::export]]
Eigen::MatrixXd thresh_C(Eigen::MatrixXd sigmaY, Eigen::MatrixXd Y, Eigen::MatrixXd Y_mean, int k, int n, int p, double deltafinal){
  double theta;
  double lambda;
  double nn=n;
  double pp=p;
  Eigen::MatrixXd lam=Eigen::MatrixXd::Zero(p,p);
  for(int i=0;i<p;i++){
    for(int j=0;j<p;j++){
      theta=0;
      for(int t=0;t<n-k;t++)
        theta = theta+pow(((Y(i,t+k)-Y_mean(i,0))*(Y(j,t)-Y_mean(j,0))-sigmaY(i,j)),2);
      theta = theta/(nn);
      lambda = deltafinal*sqrt(theta*log(pp)/nn);
      lam(i,j)=lambda;
      if(abs(sigmaY(i,j))<lambda)
        sigmaY(i,j)=0;
    }
  }
  return sigmaY;
}


// [[Rcpp::export]]
SEXP MatMult(Eigen::MatrixXd A, Eigen::MatrixXd B){
  Eigen::MatrixXd C = A * B;
  
  return Rcpp::wrap(C);
}