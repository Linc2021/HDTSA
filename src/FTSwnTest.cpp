#include <RcppEigen.h>
#include <Rcpp.h>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <random>
#include <ctime>
#include <cmath>
#include "testtools.h"

// [[Rcpp::depends(RcppEigen,Rcpp)]]
using namespace std;
using namespace Eigen;
using namespace Rcpp;

 
// // [[Rcpp::export]]
// Eigen::MatrixXd WN_XiC(int n, int k, int p, int B, double bn, int ken_sign) {
//   Eigen::MatrixXd kenel = Eigen::MatrixXd::Ones(n - k, n - k);
//   
//   // Construct Kernel Matrix based on selected type
//   if (ken_sign == 1) { // Quadratic Spectral (QS) kernel
//     for (int i = 0; i < n - k; i++) {
//       for (int j = 0; j < n - k; j++) {
//         if (i != j) {
//           double temp = double(i - j) / bn;
//           kenel(i, j) = 25 / (12 * M_PI * M_PI * temp * temp) * (sin(6 * M_PI * temp / 5) / (6 * M_PI * temp / 5) - cos(6 * M_PI * temp / 5));
//         }
//       }
//     }
//   }
//   if (ken_sign == 2) { // Parzen kernel
//     for (int i = 0; i < n - k; i++) {
//       for (int j = 0; j < n - k; j++) {
//         double temp = abs(double(i - j) / bn);
//         if (temp <= 0.5)
//           kenel(i, j) = 1 - 6 * temp * temp + 6 * temp * temp * temp;
//         else if (temp <= 1)
//           kenel(i, j) = 2 * pow((1 - temp), 3);
//         else
//           kenel(i, j) = 0;
//       }
//     }
//   }
//   if (ken_sign == 3) { // Bartlett kernel
//     for (int i = 0; i < n - k; i++) {
//       for (int j = 0; j < n - k; j++) {
//         double temp = abs(double(i - j) / bn);
//         if (temp <= 1)
//           kenel(i, j) = 1 - temp;
//         else
//           kenel(i, j) = 0;
//       }
//     }
//   }
//   
//   // Generate standard normal random variables
//   static default_random_engine e(time(0));
//   static normal_distribution<double> normal(0.0, 1.0);
//   Eigen::MatrixXd Xi_temp = Eigen::MatrixXd::Zero(B, n - k);
//   for (int i = 0; i < B; i++) {
//     for (int j = 0; j < n - k; j++) {
//       Xi_temp(i, j) = normal(e);
//     }
//   }
//   
//   // Eigenvalue decomposition to impose temporal structure
//   EigenSolver<MatrixXd> eig(kenel);
//   Eigen::VectorXd EigenValue = eig.eigenvalues().real().array();
//   // Handling non-positive eigenvalues for numerical stability
//   for (int i = 0; i < n - k; i++) {
//     if (EigenValue(i) < 0) EigenValue(i) = log(double(p)) / double(n);
//   }
//   EigenValue = EigenValue.array().sqrt();
//   Eigen::MatrixXd D = EigenValue.asDiagonal();
//   Eigen::MatrixXd EigenVector = eig.eigenvectors().real();
//   
//   // Transform weights: Xi = (Q * D^0.5 * Q^T * Z)^T
//   Eigen::MatrixXd Xi = (EigenVector * D * EigenVector.transpose() * Xi_temp.transpose()).transpose();
//   return Xi;
// }

// ------------------------------ Matrix Version -----------------------

/**
 * Calculate the test statistic Tn using matrix operations.
 * Computes max absolute cross-covariance across all lags l up to L.
 */
//[[Rcpp::export]]
double TnC_mat(int n, int L, Eigen::MatrixXd Vep) {
  Eigen::VectorXd re_L1 = Eigen::VectorXd::Zero(L);
  for (int l = 1; l < L + 1; l++) {
    // Calculate Sigma at lag l: (Vep_lag0)^T * (Vep_lag_l) / (n-l)
    Eigen::MatrixXd hat_Sigma_l_i_j = Vep.topRows(n - l).transpose() * Vep.bottomRows(n - l) / double(n - l);
    double hat_Sigma_l_i_j_abs = hat_Sigma_l_i_j.array().abs().maxCoeff();
    re_L1(l - 1) = hat_Sigma_l_i_j_abs;
  }
  double tn = re_L1.maxCoeff() * sqrt(n);
  return tn;
}

/**
 * Generate the intermediate 'ft' matrix (Kronecker products of residual vectors).
 * This represents the series whose covariance structure we analyze for bandwidth selection.
 */
// [[Rcpp::export]]
Eigen::MatrixXd WN_ftC_mat(int n, int L, int p, int M, Eigen::MatrixXd Vep) {
  int dim = p * p * M * M;
  Eigen::MatrixXd ft = Eigen::MatrixXd::Zero(L * dim, n - L);
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < (n - L); j++) {
      // Block row represents interaction between time j and time j+i+1
      ft.block(i * dim, j, dim, 1) = kroneckerProduct(Vep.row(j + i + 1).transpose(), Vep.row(j).transpose());
    }
  }
  // Demean the rows
  ft = ft - (ft.rowwise().sum() / double(n - L)).replicate(1, n - L);
  return ft;
}

/**
 * Bandwidth selection for the kernel using the Matrix 'ft' as input.
 * Estimates AR(1) parameters for each series in ft to find optimal smoothing parameter bn.
 */
// [[Rcpp::export]]
double WN_bandwith_mat(Eigen::MatrixXd ft, int n, int L, int p, int M, int ken_type) {
  int kpp = L * p * p * M * M;
  double a_hat, bw = 0;
  double atemp1, atemp2;
  Eigen::VectorXd rho = Eigen::VectorXd::Zero(kpp);
  Eigen::VectorXd sig = Eigen::VectorXd::Zero(kpp);
  
  for (int i = 0; i < kpp; i++) {
    Eigen::MatrixXd ft_tmp = ft.row(i);
    // Estimate AR(1) coefficient rho
    Eigen::MatrixXd rho_tmp = (ft_tmp.leftCols(n - L - 1) * ft_tmp.rightCols(n - L - 1).transpose()).array() / 
      ((ft_tmp.leftCols(n - L - 1) * ft_tmp.leftCols(n - L - 1).transpose()).array() + 0.000001);
    rho(i) = rho_tmp(0, 0);
    // Estimate residual variance sigma^2
    Eigen::MatrixXd sig_tmp = ((ft_tmp.rightCols(n - L - 1) - rho(i) * ft_tmp.leftCols(n - L - 1)).array().square().rowwise().sum()) / double(n - L - 1);
    sig(i) = sig_tmp(0, 0);
  }
  
  atemp1 = 0.0; atemp2 = 0.0;
  if (ken_type == 1 || ken_type == 2) { // QS or Parzen
    for (int i = 0; i < kpp; i++) {
      atemp1 += 4 * rho(i) * rho(i) * pow(sig(i), 2) * pow((1 - rho(i)), -8);
      atemp2 += pow(sig(i), 2) * pow((1 - rho(i)), -4);
    }
    a_hat = atemp1 / atemp2;
    if (ken_type == 1) bw = 1.3221 * pow((n * a_hat), 0.2);
    else bw = 2.6614 * pow((n * a_hat), 0.2);
  } else if (ken_type == 3) { // Bartlett
    for (int i = 0; i < kpp; i++) {
      atemp1 += 4 * rho(i) * rho(i) * pow(sig(i), 2) * pow((1 - rho(i)), -6) * pow((1 + rho(i)), -2);
      atemp2 += pow(sig(i), 2) * pow((1 - rho(i)), -4);
    }
    a_hat = atemp1 / atemp2;
    bw = 1.1447 * pow((n * a_hat), 0.33333333);
  }
  return bw;
}

/**
 * Bootstrap the test statistic using Parametric Bootstrap (Matrix version).
 * Perturbs the original series with random weights Xi.
 */
// [[Rcpp::export]]
Eigen::MatrixXd newbootC_kernel_mat(const int n, const int L, const int p, const int M, const int B,
                             double bn, int method, Eigen::MatrixXd Vep,
                             Eigen::MatrixXd Xi_temp) {
  // Generate weight matrix (B samples)
  Eigen::MatrixXd Xi = XiC(n, L, p, B, bn, method, Xi_temp).transpose(); // (n-L) x B
  Eigen::MatrixXd Vep_temp = Vep.topRows(n - L); // (n-L)*pM
  Eigen::MatrixXd tnstar = Eigen::MatrixXd::Zero(1, B);
  Eigen::MatrixXd re_L1 = Eigen::MatrixXd::Zero(L, B);
  
  for (int b = 0; b < B; b++) {
    Eigen::MatrixXd w_b = Xi.col(b).replicate(1, p * M); // Expand weight vector to match residual matrix dim (n-L)*pM
    Eigen::MatrixXd Vep_new = Vep_temp.array() * w_b.array(); // Multiply residuals by bootstrap weights (n-L)*pM
    
    for (int l = 1; l < L + 1; l++) {
      // Compute perturbed cross-covariance at lag l
      Eigen::MatrixXd hat_Sigma_l_i_j = Vep_new.transpose() * Vep.middleRows(l, n - L) -
        (Xi.col(b).array().sum() / double(n - L)) * Vep_temp.transpose() * Vep.middleRows(l, n - L); //pM*pM
      re_L1(l - 1, b) = hat_Sigma_l_i_j.array().abs().maxCoeff();
    }
  }
  tnstar.row(0) = re_L1.colwise().maxCoeff() / sqrt(double(n - L));
  return (tnstar);
}

// -------------------- Vectorization Version ---------------------------

/**
 * Calculate the test statistic Tn using element-wise loops.
 * More memory-efficient than the matrix version for large dimensions.
 */
//[[Rcpp::export]]
double TnC_vec(int n, int L, Eigen::MatrixXd Vep) {
  int pM = Vep.cols();
  double Tn = 0.0;
  for (int l = 1; l < L + 1; l++) {
    for (int pj = 0; pj < pM; pj++) {
      for (int pk = 0; pk < pM; pk++) {
        // Compute single covariance element for lag l between channel pj and pk
        Eigen::MatrixXd Sigma_lijuv = Vep.block(0, pk, n - l, 1).transpose() * Vep.block(l, pj, n - l, 1) / double(n - l);
        Tn = max(Tn, abs(Sigma_lijuv(0, 0)));
      } //pM
    } //pM
  } // L
  return Tn * sqrt(n);
}

/**
 * Bandwidth selection using the Vector approach.
 * Instead of taking a giant 'ft' matrix, it computes ft_tmp on the fly for each pair (pj, pk).
 */
// [[Rcpp::export]]
double WN_bandwith_vec(Eigen::MatrixXd Vep, int n, int L, int p, int M, int ken_type) {
  double a_hat, bw = 0.0;
  double atemp1 = 0.0, atemp2 = 0.0;
  
  for (int l = 1; l < L + 1; l++) {
    for (int pj = 0; pj < p * M; pj++) {
      for (int pk = 0; pk < p * M; pk++) {
        // Generate ft element-wise on the fly
        Eigen::MatrixXd ft_tmp = (Vep.block(0, pj, n - L, 1).array() * Vep.block(l, pk, n - L, 1).array()).transpose();
        ft_tmp = ft_tmp - (ft_tmp.rowwise().sum() / double(n - L)).replicate(1, n - L);
        
        // AR(1) estimation
        Eigen::MatrixXd rho_tmp = (ft_tmp.leftCols(n - L - 1) * ft_tmp.rightCols(n - L - 1).transpose()).array() / 
          ((ft_tmp.leftCols(n - L - 1) * ft_tmp.leftCols(n - L - 1).transpose()).array() + 0.000001);
        double rho = rho_tmp(0, 0);
        Eigen::MatrixXd sig_tmp = ((ft_tmp.rightCols(n - L - 1) - rho * ft_tmp.leftCols(n - L - 1)).array().square().rowwise().sum()) / double(n - L - 1);
        double sig = sig_tmp(0, 0);
        
        // Accumulate statistics based on kernel type
        if (ken_type == 1 || ken_type == 2) {
          atemp1 += 4 * rho * rho * pow(sig, 2) * pow((1 - rho), -8);
          atemp2 += pow(sig, 2) * pow((1 - rho), -4);
        } else if (ken_type == 3) {
          atemp1 += 4 * rho * rho * pow(sig, 2) * pow((1 - rho), -6) * pow((1 + rho), -2);
          atemp2 += pow(sig, 2) * pow((1 - rho), -4);
        } 
      } //pk
    } //pj
  }// l
  
  a_hat = atemp1 / atemp2;
  if (ken_type == 1) bw = 1.3221 * pow(((n - L) * a_hat), 0.2);
  else if (ken_type == 2) bw = 2.6614 * pow(((n - L) * a_hat), 0.2);
  else if (ken_type == 3) bw = 1.1447 * pow(((n - L) * a_hat), 0.33333333);
  
  return bw;
}

/**
 * Bootstrap the test statistic using the Vectorized approach.
 * Calculates bootstrap statistics without storing massive intermediate Kronecker product matrices.
 */
// [[Rcpp::export]]
Eigen::MatrixXd newbootC_kernel_vec(const int n, const int L, const int p, const int M, const int B,
                             double bn, int method, Eigen::MatrixXd Vep,
                             Eigen::MatrixXd Xi_temp) {
  Eigen::MatrixXd Xi = XiC(n, L, p, B, bn, method, Xi_temp).transpose(); // (n-L) x B
  Eigen::MatrixXd tnstar = Eigen::MatrixXd::Zero(1, B);
  
  for (int l = 1; l < L + 1; l++) {
    for (int pj = 0; pj < p * M; pj++) {
      for (int pk = 0; pk < p * M; pk++) {
        // Reconstruct ft series for lag l
        Eigen::MatrixXd ft_tmp = (Vep.block(0, pj, n - L, 1).array() * Vep.block(l, pk, n - L, 1).array()).transpose();
        ft_tmp = ft_tmp - (ft_tmp.rowwise().sum() / double(n - L)).replicate(1, n - L);
        
        // Directly multiply by weight matrix to get bootstrap statistics for all B samples
        Eigen::MatrixXd tnstar_tmp = ((ft_tmp * Xi).array() / sqrt(double(n - L))).abs();
        // Update global maximum over all lags and pairs
        tnstar = tnstar.array().max(tnstar_tmp.array());
      }
    }
  }
  return (tnstar);
}