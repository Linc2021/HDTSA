#include <RcppEigen.h>
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include "testtools.h"

// [[Rcpp::depends(RcppEigen,Rcpp)]]

namespace {

double kernel_value(const double x, const int method) {
  if (method == 1) {
    if (x == 0.0) return 1.0;
    const double z = 6.0 * M_PI * x / 5.0;
    return 25.0 / (12.0 * M_PI * M_PI * x * x) *
      (std::sin(z) / z - std::cos(z));
  }
  // The original Parzen and Bartlett branches use abs((i - j) / bn).
  const double ax = std::abs(x);
  if (method == 2) {
    if (ax <= 0.5) return 1.0 - 6.0 * ax * ax + 6.0 * ax * ax * ax;
    if (ax <= 1.0) return 2.0 * std::pow(1.0 - ax, 3.0);
    return 0.0;
  }
  if (method == 3) return ax <= 1.0 ? 1.0 - ax : 0.0;
  Rcpp::stop("Unknown kernel method.");
}
  
  Eigen::MatrixXd kernel_matrix(const Eigen::Index size, const double bandwidth,
                                const int method) {
    Eigen::MatrixXd kernel = Eigen::MatrixXd::Identity(size, size);
    for (Eigen::Index i = 0; i < size; ++i) {
      for (Eigen::Index j = i + 1; j < size; ++j) {
        const double value = kernel_value(
          (static_cast<double>(i) - static_cast<double>(j)) / bandwidth,
          method
        );
        kernel(i, j) = value;
        kernel(j, i) = value;
      }
    }
    return kernel;
  }
  
  Eigen::MatrixXd correlate_normals(const Eigen::MatrixXd& normals,
                                    const double bandwidth,
                                    const int method) {
    const Eigen::MatrixXd kernel =
      kernel_matrix(normals.cols(), bandwidth, method);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(kernel);
    if (eig.info() != Eigen::Success)
      Rcpp::stop("Eigendecomposition of the multiplier kernel failed.");
    
    Eigen::VectorXd eigenvalues = eig.eigenvalues();
    for (Eigen::Index i = 0; i < eigenvalues.size(); ++i) {
      // Supported kernels are positive semidefinite. Negative values here are
      // floating-point roundoff and should be projected to zero.
      if (eigenvalues(i) < 0.0) eigenvalues(i) = 0.0;
    }
    Eigen::MatrixXd transformed = normals * eig.eigenvectors();
    transformed.array().rowwise() *= eigenvalues.cwiseSqrt().transpose().array();
    return transformed * eig.eigenvectors().transpose();
  }
  
  class BandwidthAccumulator {
  public:
    BandwidthAccumulator(const int method, const Eigen::Index n_tilde)
      : method_(method), n_tilde_(n_tilde), numerator_(0.0L), denominator_(0.0L) {}
    
    void add(const double rho, const double innovation_variance) {
      const double one_minus_rho = 1.0 - rho;
      if (!std::isfinite(rho) || !std::isfinite(innovation_variance) ||
          !(innovation_variance >= 0.0) ||
          std::abs(one_minus_rho) <= std::numeric_limits<double>::epsilon()) {
        return;
      }
      const long double sigma4 =
        static_cast<long double>(innovation_variance) * innovation_variance;
      denominator_ += sigma4 * std::pow(one_minus_rho, -4.0);
      if (method_ == 3) {
        const double one_plus_rho = 1.0 + rho;
        if (std::abs(one_plus_rho) <= std::numeric_limits<double>::epsilon()) return;
        numerator_ += 4.0L * rho * rho * sigma4 *
          std::pow(one_minus_rho, -6.0) * std::pow(one_plus_rho, -2.0);
      } else {
        numerator_ += 4.0L * rho * rho * sigma4 *
          std::pow(one_minus_rho, -8.0);
      }
    }
    
    double result() const {
      if (!(denominator_ > 0.0L) || !(numerator_ >= 0.0L))
        Rcpp::stop("Unable to estimate a finite positive bandwidth from the data.");
      const double a_hat = static_cast<double>(numerator_ / denominator_);
      double bandwidth;
      if (method_ == 1) {
        bandwidth = 1.3221 * std::pow(static_cast<double>(n_tilde_) * a_hat, 0.2);
      } else if (method_ == 2) {
        bandwidth = 2.6614 * std::pow(static_cast<double>(n_tilde_) * a_hat, 0.2);
      } else if (method_ == 3) {
        bandwidth = 1.1447 *
          std::pow(static_cast<double>(n_tilde_) * a_hat, 1.0 / 3.0);
      } else {
        Rcpp::stop("Unknown kernel method.");
      }
      if (!(bandwidth > 0.0) || !std::isfinite(bandwidth))
        Rcpp::stop("The estimated bandwidth is not finite and positive.");
      return bandwidth;
    }
    
  private:
    int method_;
    Eigen::Index n_tilde_;
    long double numerator_;
    long double denominator_;
  };
  
} // namespace


Eigen::MatrixXd XiC(int n, int k, int p, int B, double bn, int ken_sign,
                    Eigen::MatrixXd Xi_temp) {
  (void) n;
  (void) k;
  (void) p;
  (void) B;
  return correlate_normals(Xi_temp, bn, ken_sign);
}


double bandwith(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Xj,
                const int K, const int p, const int d, const int method) {
  const Eigen::Index n_tilde = X.rows() - K;
  if (X.cols() != p) Rcpp::stop("p does not match ncol(X).");
  if (Xj.rows() != X.rows() || Xj.cols() != d)
    Rcpp::stop("Xj dimensions do not match nrow(X) and d.");
  BandwidthAccumulator accumulator(method, n_tilde);
  
  for (int lag = 1; lag <= K; ++lag) {
    const Eigen::MatrixXd future = X.middleRows(lag, n_tilde);
    const Eigen::MatrixXd past = Xj.topRows(n_tilde);
    const Eigen::MatrixXd future_left = future.topRows(n_tilde - 1);
    const Eigen::MatrixXd future_right = future.bottomRows(n_tilde - 1);
    const Eigen::MatrixXd past_left = past.topRows(n_tilde - 1);
    const Eigen::MatrixXd past_right = past.bottomRows(n_tilde - 1);
    
    const Eigen::MatrixXd means =
      (future.transpose() * past) / static_cast<double>(n_tilde);
    const Eigen::MatrixXd sum_left = future_left.transpose() * past_left;
    const Eigen::MatrixXd sum_right = future_right.transpose() * past_right;
    const Eigen::MatrixXd square_left =
      future_left.array().square().matrix().transpose() *
      past_left.array().square().matrix();
    const Eigen::MatrixXd square_right =
      future_right.array().square().matrix().transpose() *
      past_right.array().square().matrix();
    const Eigen::MatrixXd raw_cross =
      (future_left.array() * future_right.array()).matrix().transpose() *
      (past_left.array() * past_right.array()).matrix();
    const double count = static_cast<double>(n_tilde - 1);
    const Eigen::MatrixXd centered_left_square =
      square_left.array() - 2.0 * means.array() * sum_left.array() +
      count * means.array().square();
    const Eigen::MatrixXd centered_right_square =
      square_right.array() - 2.0 * means.array() * sum_right.array() +
      count * means.array().square();
    const Eigen::MatrixXd centered_cross =
      raw_cross.array() - means.array() * sum_left.array() -
      means.array() * sum_right.array() + count * means.array().square();
    
    for (Eigen::Index col = 0; col < d; ++col) {
      for (Eigen::Index row = 0; row < p; ++row) {
        const double ar_denom = centered_left_square(row, col);
        if (!(ar_denom > std::numeric_limits<double>::epsilon())) continue;
        const double cross = centered_cross(row, col);
        const double rho = cross / ar_denom;
        const double innovation_variance = std::max(
          0.0,
          (centered_right_square(row, col) - 2.0 * rho * cross +
            rho * rho * ar_denom) / count
        );
        accumulator.add(rho, innovation_variance);
      }
    }
  }
  return accumulator.result();
}
