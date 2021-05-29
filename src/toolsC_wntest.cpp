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
Eigen::MatrixXd WN_XiC(int n, int k,int p, int B, double bn, int ken_sign){
  Eigen::MatrixXd kenel = Eigen::MatrixXd::Ones(n-k,n-k);
  if(ken_sign==1){
    for(int i=0; i<n-k; i++){
      for(int j=0; j<n-k; j++){
        if(i!=j){
          double temp =double(i-j)/bn;
          kenel(i,j) = 25/(12*PI*PI*temp*temp) * (sin(6*PI*temp/5)/(6*PI*temp/5) - cos(6*PI*temp/5));
        }
      }
    }
  }
  if(ken_sign==2){
    for(int i=0; i<n-k; i++){
      for(int j=0; j<n-k; j++){
        double temp =abs(double(i-j)/bn);
        if(temp<=0.5)
          kenel(i,j) = 1-6*temp*temp+6*temp*temp*temp;
        else if(temp<=1)
          kenel(i,j) = 2*pow((1-temp),3);
        else
          kenel(i,j) = 0;
      }
    }
  }
  if(ken_sign==3){
    for(int i=0; i<n-k; i++){
      for(int j=0; j<n-k; j++){
        double temp =abs(double(i-j)/bn);
        if(temp<=1)
          kenel(i,j) = 1-temp;
        else
          kenel(i,j) = 0;
      }
    }
  }
  static default_random_engine e(time(0));
  static normal_distribution<double> normal(0.0,1.0);
  Eigen::MatrixXd Xi_temp=Eigen::MatrixXd::Zero(B,n-k);
  for(int i=0; i<B; i++){
    for(int j=0; j<n-k; j++){
      Xi_temp(i,j) = normal(e);    
    }
  }
  EigenSolver<MatrixXd> eig(kenel);
  Eigen::VectorXd EigenValue =  eig.eigenvalues().real().array();
  for(int i=0;i<n-k;i++){
    if(EigenValue(i)<0)EigenValue(i)=log(double(p))/double(n);
  }
  EigenValue = EigenValue.array().sqrt();
  Eigen::MatrixXd D = EigenValue.asDiagonal();
  Eigen::MatrixXd EigenVector = eig.eigenvectors().real();
  Eigen::MatrixXd Xi = (EigenVector*D*EigenVector.transpose()*Xi_temp.transpose()).transpose();
  return Xi;
}

// [[Rcpp::export]]
double WN_bandwith(Eigen::MatrixXd ft, int k,int p, int ken_type){
  int n = ft.cols();
  int kpp = k*p*p;
  double a_hat,bw=0;
  double atemp1,atemp2;
  Eigen::VectorXd rho = (ft.leftCols(n-1)*ft.rightCols(n-1).transpose()).diagonal().array()/(ft.leftCols(n-1)*ft.leftCols(n-1).transpose()).diagonal().array();
  Eigen::VectorXd sig = ((ft.rightCols(n-1)-rho.asDiagonal()*ft.leftCols(n-1)).array().square().rowwise().sum())/double(n-1);
  //VectorXd bw_list = VectorXd::Zero(kpd*(kpd-1)/2);
  atemp1 = 0.0;
  atemp2 = atemp1;
  if(ken_type==1||ken_type==2){
    for(int i=0;i<kpp;i++){
      atemp1 += 4*rho(i)*rho(i)*pow(sig(i),2)*pow((1-rho(i)),-8);
      atemp2 += pow(sig(i),2)*pow((1-rho(i)),-4);
    }
    a_hat = atemp1/atemp2;
    if(ken_type==1) bw = 1.3221*pow(((n-k)*a_hat),0.2);
    else bw = 2.6614*pow(((n-k)*a_hat),0.2);
  }
  else if(ken_type==3){
    for(int i=0;i<kpp;i++){
      atemp1 += 4*rho(i)*rho(i)*pow(sig(i),2)*pow((1-rho(i)),-6)*pow((1+rho(i)),-2);
      atemp2 += pow(sig(i),2)*pow((1-rho(i)),-4);
    }
    a_hat = atemp1/atemp2;
    bw = 1.1447*pow(((n-k)*a_hat),0.33333333);
  }
  return bw;
}



// [[Rcpp::export]]
std::vector<double> WN_bootc(const int n, const int k, const int p, const int B,  
                                 double bn,int method,Eigen::MatrixXd ft, Eigen::MatrixXd X, Eigen::VectorXd sigma_zero){
  
  // random samples follows a multivariate gaussian distribution N(0, Jn):B*kpd
  Eigen::MatrixXd Xi= WN_XiC(n, k, p, B, bn, method);  //Xi(B,n-k)
  Eigen::VectorXd Wk = kroneckerProduct(sigma_zero,sigma_zero);
  Eigen::VectorXd Ik = VectorXd::Ones(k);
  Eigen::VectorXd W = kroneckerProduct(Ik,Wk);
  //Rcout<<2<<"\n";
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


// std::vector<double> newbootc_wn(const int n, const int k, const int p, const int B,  
//                                 const int M, MatrixXd X, VectorXd sigma_zero){
//   
//   int Q = floor((n-k)/M);
//   VectorXd Wk = kroneckerProduct(sigma_zero,sigma_zero);
//   //Rcout<<2<<"\n";
//   VectorXd Ik = VectorXd::Ones(k);
//   //Rcout<<2<<"\n";
//   VectorXd W = kroneckerProduct(Ik,Wk);
//   //Rcout<<2<<"\n";
//   MatrixXd ft = ftC_wn(n, k, p, Q, M, X);
//   //Rcout<<2<<"\n";
//   for(int i=0; i<(Q*M); i++){
//     ft.col(i) = ft.col(i).array()*W.array();
//   }
//   //Rcout<<2<<"\n";
//   // random samples follows a multivariate gaussian distribution N(0, Jn):B*kpd
//   MatrixXd eq = eqC(Q,B,M);
//   MatrixXd samples = eq * ft.leftCols(Q*M).transpose()/ sqrt(Q*M);
//   //Rcout<<2<<"\n";
//   //Rcout<<2<<"\n";
//   
//   std::vector<double> GnStar(B);
//   for(int b=0; b<B; b++){
//     
//     GnStar[b] = samples.row(b).array().abs().maxCoeff();
//     
//   }
//   sort(GnStar.begin(),GnStar.end());
//   return GnStar;
// }

// // [[Rcpp::export]]
// Rcpp::IntegerMatrix leaveoneC_wn(int n, int k, int p, int B, int M, int cv_num, double alpha,
//                                  Rcpp::IntegerVector sample, MatrixXd X){
//   
//   int loc;
//   double Tn;
//   List Tn_tools;
//   Rcpp::IntegerMatrix reject(cv_num, M);
//   MatrixXd X_tmp;
//   std::vector<double> Gnstar;
//   for(int i=0 ;i<cv_num; i++){
//     loc = sample(i);
//     X_tmp = oneout(n, p, loc, X);
//     Tn_tools = teststatC _wn(X_tmp, n-1, p, k);
//     Tn = Tn_tools["Tn"];
//     //if(i%10==0){Rcout<<i;}
//     VectorXd sigma_zero = Tn_tools["sigmaO"];
//     for(int m=0; m<M; m++){
//       Gnstar = newbootc_wn(n-1, k, p, B, m+1, X_tmp, sigma_zero);
//       reject(i, m) = (Tn>Gnstar[floor(B*(1-alpha))-1])? 1:0;
//       
//     }
//   }
//   return reject;
// }  
