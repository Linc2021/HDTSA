#ifndef HDTSA_TESTTOOLS_H
#define HDTSA_TESTTOOLS_H

#include <RcppEigen.h>

Eigen::MatrixXd XiC(int n, int k, int p, int B, double bn, int ken_sign,
                    Eigen::MatrixXd Xi_temp);

// Common low-memory bandwidth estimator. For WN use Xj = X and d = p; for the
// martingale-difference test Xj is the transformed predictor matrix.
double bandwith(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Xj,
                int k, int p, int d, int ken_type);

#endif
