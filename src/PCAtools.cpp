
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
Eigen::MatrixXd thresh_C(Eigen::MatrixXd mat, double delta){
  //double threshold = lambda * sqrt(log(p) / n);
  for (int i = 0; i < mat.rows(); i++) {
    for (int j = 0; j < mat.cols(); j++) {
      if (mat(i, j) < delta) {
        mat(i, j) = 0;
      }
    }
  }
  return mat;
}


// [[Rcpp::export]]
SEXP MatMult(Eigen::MatrixXd A, Eigen::MatrixXd B){
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}