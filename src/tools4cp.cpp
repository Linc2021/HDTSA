#include <RcppEigen.h>


// [[Rcpp::depends(RcppEigen)]]

using namespace Eigen;


// [[Rcpp::export]]
Eigen::VectorXd minor_P(Eigen::MatrixXd Wr, Eigen::MatrixXd Ws, int d1, int d2){

  Eigen::VectorXd P(d1*d1*d2*d2);
  int sss = 0;
  for(int l = 0; l < d2;l++){
    for(int k = 0; k < d2;k++){
      for(int j = 0; j < d1;j++){
        for(int i = 0; i < d1;i++){
          P(sss) = Wr(i,k) * Ws(j,l) + Ws(i,k) * Wr(j,l) - Wr(i,l) * Ws(j,k) - Ws(i,l) * Wr(j,k);
          sss++;
        }
      }
    }
  }

  return P;

}

// [[Rcpp::export]]
Eigen::MatrixXd Vech2Mat_new(Eigen::VectorXd P, int d){

  Eigen::MatrixXd M  = Eigen::MatrixXd::Zero(d,d);

  int k = 0;
  for(int j = 0; j < d;j++){
    for(int i = j; i < d;i++){
      M(i,j) = P(k)/2;
      k++;
    }
  }

  M = M + M.transpose();

  return M;

}
