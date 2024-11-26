#include <RcppEigen.h>
#include <Rcpp.h>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <random>
#include <ctime>
// [[Rcpp::depends(RcppEigen,Rcpp)]]
using namespace std;
using namespace Eigen;
using namespace Rcpp;

Eigen::MatrixXd XiC(int n, int k,int p, int B, double bn, int ken_sign, Eigen::MatrixXd Xi_temp);

double bandwith(Eigen::MatrixXd ft, int k,int p, int d, int ken_type);



