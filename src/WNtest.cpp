#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
// [[Rcpp::depends(Rcpp)]]
#include <iostream>
#include <algorithm>
#include <math.h>
#include <random>
#include <ctime>
#include "testtools.h"
using namespace std;
using namespace Eigen;
using namespace Rcpp;



// [[Rcpp::export]]
Rcpp::List WN_teststatC(Eigen::MatrixXd X, int n, int p, int k){
  double nn = n;
  Eigen::VectorXd Tnk = Eigen::VectorXd::Zero(k);
  Eigen::MatrixXd X_meanVec = (X.colwise().sum()/nn);
  Eigen::MatrixXd X_mean = X_meanVec.array().replicate(n,1);
  X = X-X_mean;
  Eigen::MatrixXd sigma_zero = ((X.transpose() * X/nn).diagonal()).cwiseSqrt().cwiseInverse().asDiagonal();
  for (int i=1; i<k+1; i++){
    double jj = n-i;
    Eigen::MatrixXd sigma_k = ((X.bottomRows(n-i)).transpose() * X.topRows(n-i)) / jj;
    //Rcout<<1<<"\n";
    Tnk(i-1) =  (sigma_zero*sigma_k*sigma_zero*sqrt(nn)).array().abs().maxCoeff();
    //Rcout<<1<<"\n";
  }
  double Tn = Tnk.array().maxCoeff();
  //Rcout<<1<<"\n";
  return List::create(Named("Tn") = Tn, Named("sigma_zero") = sigma_zero.diagonal().array(), Named("X_mean") = X_meanVec);
}

// [[Rcpp::export]]
Eigen::MatrixXd WN_ftC(int n, int k, int p, Eigen::MatrixXd X, Eigen::MatrixXd X_mean){
  
  int dim = p*p;
  Eigen::MatrixXd ft = MatrixXd::Zero(k*p*p, n-k);
  X = X - X_mean.replicate(n,1);;
  for(int i=0; i<k; i++){
    for(int j=0; j<(n-k); j++){
      ft.block(i*dim,j,dim,1) = kroneckerProduct(X.row(j+i+1).transpose(), X.row(j).transpose());
    }
  }
  ft = ft - (ft.rowwise().sum()/double(n-k)).replicate(1,n-k);
  return ft;
}


// [[Rcpp::export]]
std::vector<double> WN_bootc(const int n, const int k, const int p, const int B,
                             double bn,int method,Eigen::MatrixXd ft,
                             Eigen::MatrixXd X, Eigen::VectorXd sigma_zero,
                             Eigen::MatrixXd Xi_temp
                                 ){
  
  // random samples follows a multivariate gaussian distribution N(0, Jn):B*kpd
  Eigen::MatrixXd Xi= XiC(n, k, p, B, bn, method, Xi_temp);  //Xi(B,n-k)
  Eigen::VectorXd Wk = kroneckerProduct(sigma_zero,sigma_zero);
  Eigen::VectorXd Ik = VectorXd::Ones(k);
  Eigen::VectorXd W = kroneckerProduct(Ik,Wk);
  for(int i=0; i<n-k; i++){
    ft.col(i) = ft.col(i).array()*W.array();
  }
  Eigen::MatrixXd samples = Xi * ft.transpose()/sqrt(double(n-k));
  std::vector<double> GnStar(B);
  for(int b=0; b<B; b++){
    
    GnStar[b] = samples.row(b).array().abs().maxCoeff();
    
  }
  sort(GnStar.begin(),GnStar.end());
  return GnStar;
}


