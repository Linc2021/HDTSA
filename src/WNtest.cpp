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
#include <limits>
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
    // Equation (2) of Chang, Yao and Zhou (2017) uses 1 / n, not
    // 1 / (n - lag), for the sample lag covariance.
    Eigen::MatrixXd sigma_k = ((X.bottomRows(n-i)).transpose() * X.topRows(n-i)) / nn;
    //Rcout<<1<<"\n";
    Tnk(i-1) =  (sigma_zero*sigma_k*sigma_zero*sqrt(nn)).array().abs().maxCoeff();
    //Rcout<<1<<"\n";
  }
  double Tn = Tnk.array().maxCoeff();
  //Rcout<<1<<"\n";
  return List::create(Named("Tn") = Tn, Named("sigma_zero") = sigma_zero.diagonal().array(), Named("X_mean") = X_meanVec);
}

namespace {

// Build only the requested coordinate block for the low-memory bootstrap.
Eigen::MatrixXd WN_ftC(const Eigen::MatrixXd& X, const int lag,
                       const Eigen::Index first, const Eigen::Index count,
                       const Eigen::Index n_tilde, const Eigen::Index p) {
  Eigen::MatrixXd out(count, n_tilde);
  if (first == 0 && count == p * p) {
    for (Eigen::Index time = 0; time < n_tilde; ++time) {
      out.col(time) = kroneckerProduct(
        X.row(time + lag).transpose(), X.row(time).transpose()
      );
    }
  } else {
    for (Eigen::Index row = 0; row < count; ++row) {
      const Eigen::Index coordinate = first + row;
      const Eigen::Index future_component = coordinate / p;
      const Eigen::Index past_component = coordinate % p;
      out.row(row) =
        (X.col(future_component).segment(lag, n_tilde).array() *
        X.col(past_component).head(n_tilde).array()).matrix().transpose();
    }
  }
  out.colwise() -= out.rowwise().mean();
  return out;
}
  
} // namespace


// Fused, low-memory implementation of equations (4), (8), (10), and (11).
// It never materializes the K p^2 by (n-K) matrix f or the B by K p^2
// bootstrap matrix. Instead, coordinates are generated and multiplied in blocks.
// [[Rcpp::export]]
Rcpp::List WN_LinfC(Eigen::MatrixXd X, const int K, const int B,
                    const int method, Eigen::MatrixXd Xi_temp,
                    int block_size = 0) {
  const Eigen::Index n = X.rows();
  const Eigen::Index p = X.cols();
  if (K < 1 || static_cast<Eigen::Index>(K) >= n - 1) {
    Rcpp::stop("lag.k must be positive and smaller than n - 1.");
  }
  if (B < 1) Rcpp::stop("B must be positive.");
  if (method < 1 || method > 3) Rcpp::stop("Unknown kernel method.");
  if (p > 0 && p > std::numeric_limits<Eigen::Index>::max() / p) {
    Rcpp::stop("p is too large for p^2 coordinate indexing.");
  }
  
  X.rowwise() -= X.colwise().mean();
  const Eigen::VectorXd variances =
    (X.array().square().colwise().sum() / static_cast<double>(n))
      .matrix().transpose();
    for (Eigen::Index i = 0; i < p; ++i) {
      if (!(variances(i) > 0.0) || !std::isfinite(variances(i))) {
        Rcpp::stop("Every component must have a finite positive sample variance.");
      }
    }
    const Eigen::VectorXd inverse_sd = variances.array().sqrt().inverse();
    
    double statistic = 0.0;
    for (int lag = 1; lag <= K; ++lag) {
      const Eigen::MatrixXd covariance =
        X.bottomRows(n - lag).transpose() * X.topRows(n - lag) /
          static_cast<double>(n);
      const Eigen::MatrixXd correlation =
        covariance.array() *
        (inverse_sd * inverse_sd.transpose()).array();
      statistic = std::max(statistic,
                           std::sqrt(static_cast<double>(n)) *
                             correlation.cwiseAbs().maxCoeff());
    }
    
    const Eigen::Index n_tilde = n - K;
    if (Xi_temp.rows() != B || Xi_temp.cols() != n_tilde)
      Rcpp::stop("boot_nomal must be a B by (n - lag.k) matrix.");
    const Eigen::Index p2 = p * p;
    if (block_size <= 0) {
      const double bytes_per_coordinate =
        8.0 * (static_cast<double>(B) + static_cast<double>(n_tilde));
      const Eigen::Index memory_limited = static_cast<Eigen::Index>(
        std::max(1.0, std::floor((128.0 * 1024.0 * 1024.0) / bytes_per_coordinate))
      );
      // Use as much of a bounded 128 MiB workspace as possible. Larger blocks are
      // materially faster because they turn many small GEMMs into a few BLAS-3
      // operations, while memory remains independent of p^2.
      block_size = static_cast<int>(std::min<Eigen::Index>(p2, memory_limited));
    }
    const Eigen::Index block = std::max<Eigen::Index>(1, block_size);
    const double bandwidth = bandwith(X, X, K, p, p, method);
    Eigen::MatrixXd multipliers =
      XiC(n, K, p, B, bandwidth, method, Xi_temp);
    Eigen::VectorXd maxima = Eigen::VectorXd::Zero(B);
    
    for (int lag = 1; lag <= K; ++lag) {
      for (Eigen::Index first = 0; first < p2; first += block) {
        const Eigen::Index count = std::min(block, p2 - first);
        Eigen::MatrixXd f = WN_ftC(X, lag, first, count, n_tilde, p);
        Eigen::VectorXd coordinate_scale(count);
        for (Eigen::Index row = 0; row < count; ++row) {
          const Eigen::Index coordinate = first + row;
          coordinate_scale(row) =
            inverse_sd(coordinate / p) * inverse_sd(coordinate % p) /
              std::sqrt(static_cast<double>(n_tilde));
        }
        f.array().colwise() *= coordinate_scale.array();
        Eigen::MatrixXd samples(B, count);
        samples.noalias() = multipliers * f.transpose();
        for (int bootstrap = 0; bootstrap < B; ++bootstrap) {
          maxima(bootstrap) = std::max(
            maxima(bootstrap), samples.row(bootstrap).cwiseAbs().maxCoeff()
          );
        }
      }
    }
    
    std::vector<double> sorted_maxima(maxima.data(), maxima.data() + maxima.size());
    std::sort(sorted_maxima.begin(), sorted_maxima.end());
    
    return Rcpp::List::create(
      Rcpp::Named("Tn") = statistic,
      Rcpp::Named("Gnstar") = sorted_maxima,
      Rcpp::Named("bandwidth") = bandwidth,
      Rcpp::Named("block.size") = block
    );
}

// [[Rcpp::export]]
arma::mat resampling(arma::mat X, int n, int p, int B, int tau) {
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
      arma::uvec indices = arma::join_cols(indices1, indices2);
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
  return arma::sort(Hn_B, "ascend", 1);
}
