#include <RcppEigen.h>
#include <Rcpp.h>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <random>
#include <ctime>
#include "testtools.h"
// [[Rcpp::depends(RcppEigen,Rcpp)]]
using namespace std;
using namespace Eigen;
using namespace Rcpp;


// [[Rcpp::export]]
double MartG_TestStatC(int n, int k,  Eigen::MatrixXd X, Eigen::MatrixXd Xj){
  
  
  Eigen::VectorXd GammaMax = Eigen::VectorXd::Zero(k);
  for(int j=1; j<k+1; j++){
    double jj=j;
    Eigen::MatrixXd GammajHat = (Xj.topRows(n-j)).transpose() * X.bottomRows(n-j) / (n-jj);
    GammaMax(j-1) = GammajHat.array().square().maxCoeff();
  }
  
  return n * GammaMax.sum();
  
}

// [[Rcpp::export]]
Eigen::MatrixXd MartG_ftC(int n, int k, int p, int d, Eigen::MatrixXd X, Eigen::MatrixXd Xj){
  
  int dim = p*d;
  Eigen::MatrixXd ft = Eigen::MatrixXd::Zero(k*p*d, n-k);
  for(int i=0; i<k; i++){
    for(int j=0; j<(n-k); j++){
      ft.block(i*dim,j,dim,1) = kroneckerProduct(Xj.row(j).transpose(), X.row(j+i+1).transpose());
    }
  }
  ft = ft - (ft.rowwise().sum()/(n-k)).replicate(1,n-k);
  return ft;
}


// [[Rcpp::export]]
std::vector<double> MartG_bootc(const int n, const int k, const int p, const int d, const int B,
                  double bn, int method, Eigen::MatrixXd ft, Eigen::MatrixXd Xi_temp){
  
  int dim = p*d;
  // random samples follows a multivariate gaussian distribution N(0, Jn):B*kpd
  Eigen::MatrixXd Xi= XiC(n, k, p, B, bn, method, Xi_temp);  //Xi(B,n-k)
  Eigen::MatrixXd samples = Xi * ft.transpose()/sqrt(double(n-k));
  std::vector<double> GnStar(B);
  for(int b=0; b<B; b++){
    Eigen::VectorXd GammaStar = Eigen::VectorXd::Zero(k) ;
    for(int i=0; i<k; i++){
      GammaStar(i) = samples.block(b,i*dim,1,dim).array().square().maxCoeff();
    }
    
    GnStar[b] = GammaStar.sum();
  }
  sort(GnStar.begin(),GnStar.end());
  return GnStar;
}  



