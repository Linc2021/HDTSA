
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
Eigen::MatrixXd vol_wy(Eigen::MatrixXd Y, Eigen::VectorXd Y_mean, int k0 ,int n, int p){
  Eigen::MatrixXd wy_hat=Eigen::MatrixXd::Zero(p, p);
  Eigen::VectorXd yl,yt,yt_k;
  int s;
  double yt_k_norm, yl_norm;
  for (int l=0; l<n; l++){
    yl = Y.col(l).array() - Y_mean.array();
    for (int k=1; k<=k0; k++){
      Eigen::MatrixXd sigmaY=Eigen::MatrixXd::Zero(p, p);
      for (int t=k; t<n; t++){
        s = t-k;
        yt_k = Y.col(s).array() - Y_mean.array();
        yt_k_norm = yt_k.array().square().sum();
        yl_norm = yl.array().square().sum();
        if (yt_k_norm<yl_norm){
          yt = Y.col(t).array() - Y_mean.array();
          sigmaY = sigmaY + yt * yt.transpose();
        }
      }
      sigmaY = sigmaY/double(n-k);
      sigmaY = sigmaY * sigmaY;
      wy_hat = wy_hat + sigmaY;
    }
  }
  return wy_hat;
}

// [[Rcpp::export]]
SEXP MatMult(Eigen::MatrixXd A, Eigen::MatrixXd B){
  Eigen::MatrixXd C = A * B;
  
  return Rcpp::wrap(C);
}