#include <RcppEigen.h>
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include "testtools.h"

// [[Rcpp::depends(RcppEigen,Rcpp)]]
using namespace std;
using namespace Eigen;
using namespace Rcpp;


// [[Rcpp::export]]
double MartG_TestStatC(int n, int k, Eigen::MatrixXd X,
                       Eigen::MatrixXd Xj) {
  Eigen::VectorXd GammaMax = Eigen::VectorXd::Zero(k);
  for (int j = 1; j < k + 1; ++j) {
    const double jj = j;
    const Eigen::MatrixXd GammajHat =
      Xj.topRows(n - j).transpose() * X.bottomRows(n - j) / (n - jj);
    GammaMax(j - 1) = GammajHat.array().square().maxCoeff();
  }
  return n * GammaMax.sum();
}


namespace {

// Build only the requested coordinate block for the low-memory bootstrap.
Eigen::MatrixXd MartG_ftC(const Eigen::MatrixXd& X,
                          const Eigen::MatrixXd& Xj,
                          const int lag, const Eigen::Index first,
                          const Eigen::Index count,
                          const Eigen::Index n_tilde,
                          const Eigen::Index p,
                          const Eigen::Index d) {
  Eigen::MatrixXd out(count, n_tilde);
  if (first == 0 && count == p * d) {
    for (Eigen::Index time = 0; time < n_tilde; ++time) {
      out.col(time) = kroneckerProduct(
        Xj.row(time).transpose(), X.row(time + lag).transpose()
      );
    }
  } else {
    for (Eigen::Index row = 0; row < count; ++row) {
      const Eigen::Index coordinate = first + row;
      const Eigen::Index past_component = coordinate / p;
      const Eigen::Index future_component = coordinate % p;
      out.row(row) =
        (Xj.col(past_component).head(n_tilde).array() *
        X.col(future_component).segment(lag, n_tilde).array())
         .matrix().transpose();
    }
  }
  out.colwise() -= out.rowwise().mean();
  return out;
}
  
} // namespace


// [[Rcpp::export]]
std::vector<double> MartG_bootc(const int n, const int k, const int p,
                                const int d, const int B, const int method,
                                Eigen::MatrixXd X, Eigen::MatrixXd Xj,
                                Eigen::MatrixXd Xi_temp,
                                int block_size = 0) {
  const Eigen::Index n_tilde = n - k;
  const Eigen::Index coordinates = static_cast<Eigen::Index>(p) * d;
  if (X.rows() != n || X.cols() != p || Xj.rows() != n || Xj.cols() != d)
    Rcpp::stop("X/Xj dimensions do not match n, p, and d.");
  if (k < 1 || k >= n - 1) Rcpp::stop("k must be positive and smaller than n - 1.");
  if (B < 1) Rcpp::stop("B must be positive.");
  if (Xi_temp.rows() != B || Xi_temp.cols() != n_tilde)
    Rcpp::stop("boot_nomal must be a B by (n - lag.k) matrix.");
  
  const double bn = bandwith(X, Xj, k, p, d, method);
  const Eigen::MatrixXd Xi = XiC(n, k, p, B, bn, method, Xi_temp);
  
  if (block_size <= 0) {
    const double bytes_per_coordinate =
      8.0 * (static_cast<double>(B) + static_cast<double>(n_tilde));
    const Eigen::Index memory_limited = static_cast<Eigen::Index>(
      std::max(1.0, std::floor((128.0 * 1024.0 * 1024.0) /
        bytes_per_coordinate))
    );
    block_size = static_cast<int>(
      std::min<Eigen::Index>(coordinates, memory_limited)
    );
  }
  const Eigen::Index block = std::max<Eigen::Index>(1, block_size);
  Eigen::VectorXd GnStar = Eigen::VectorXd::Zero(B);
  
  for (int lag = 1; lag <= k; ++lag) {
    Eigen::VectorXd lag_max = Eigen::VectorXd::Zero(B);
    for (Eigen::Index first = 0; first < coordinates; first += block) {
      const Eigen::Index count = std::min(block, coordinates - first);
      const Eigen::MatrixXd f =
        MartG_ftC(X, Xj, lag, first, count, n_tilde, p, d);
      const Eigen::MatrixXd samples = Xi * f.transpose() /
        std::sqrt(static_cast<double>(n_tilde));
      for (int bootstrap = 0; bootstrap < B; ++bootstrap) {
        lag_max(bootstrap) = std::max(
          lag_max(bootstrap), samples.row(bootstrap).cwiseAbs().maxCoeff()
        );
      }
    }
    GnStar.array() += lag_max.array().square();
  }
  
  std::vector<double> result(GnStar.data(), GnStar.data() + GnStar.size());
  std::sort(result.begin(), result.end());
  return result;
}
