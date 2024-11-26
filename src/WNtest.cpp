#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
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


// [[Rcpp::export]]
arma::mat resampling(arma::mat X, int n, int p, int B, int tau) {
  // You can include R code blocks in C++ files processed with sourceCpp
  // (useful for testing and development). The R code will be automatically 
  // run after the compilation.
  double constant = sqrt(2)/double(n-1);
  double Gn_1e, sigma_n1e;
  arma::mat Y, tmp;
  arma::mat Hn_B(tau, B);
  arma::vec Hne(tau, arma::fill::zeros);
  arma::vec rand_unif(n, arma::fill::zeros);
  
  for(int i = 0; i < B; i++){
    rand_unif.randu();
    arma::vec et(n, arma::fill::value(-1.0));
    et.elem(find(rand_unif > 0.5)) += 2;
    //X.print();
    //et.print();
    Y = et % X.each_col();
    tmp = Y * Y.t();
    tmp.diag().zeros();
    sigma_n1e = constant * accu(pow(tmp, 2));
    Hne.zeros();
    for(int lag=1; lag<=tau; lag++){
      arma::uvec indices1 = arma::linspace<arma::uvec>(n - lag, n-1, lag);
      arma::uvec indices2 = arma::linspace<arma::uvec>(0, n - 1- lag, n - lag);
      arma::uvec indices = join_cols(indices1, indices2);
      // if(i==0)
      //   indices.print();
      // Gn_1e = accu(tmp.submat(0 , 0, n-1-lag, n-1-lag ) % tmp.submat(lag , lag, n-1, n-1 ));
      Gn_1e = accu(tmp % tmp.submat(indices , indices));
      // Rcpp::Rcout << Gn_1e <<'\n'<< std::endl;
      if(lag > 1){
        Hne(lag-1) = Hne(lag-2) + Gn_1e/sigma_n1e;
      }
      else{
        Hne(lag-1) = Gn_1e / sigma_n1e;
      }
      
      // Hne(lag-1) = sum(Hne) + Gn_1e / sigma_n1e;
    }
    Hn_B.col(i) = Hne;
  }
  return sort(Hn_B, "ascend", 1);
}