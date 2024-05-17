// [[Rcpp::depends(RcppEigen,Rcpp)]]
#include <RcppEigen.h>
#include <Rcpp.h>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <random>
#include <ctime>


using namespace std;
using namespace Eigen;
using namespace Rcpp;

// ===================================================================
//  Kernel function
// ===================================================================

// [[Rcpp::export]]
double TaperQsC(double x)
{
  double z = 6*M_PI*abs(x)/5;
  double value;
  if(z==0.0){
    value = 1;
  }else{
    value = 3*(sin(z) - z*cos(z))/pow(z,3);
  }
  return value;
}

// [[Rcpp::export]]
double TaperBartC(double x)
{
  double z = 1-abs(x);
  
  return max(z, 0.0);
  
}

// [[Rcpp::export]]
double TaperFlatC(double x, double c)
{
  
  double z = (1-abs(x)) / (1-c);
  
  return min(max(z,0.0), 1.0);
  
}


// ===================================================================
//  compute the optimal bandwidth for the estimation of spectral density
// ===================================================================

// [[Rcpp::export]]
Rcpp::List CmpGammaC(Eigen::MatrixXd Vt)
{
  // A function to compute Autocovariance matrix sequence Gamma
  // Vt: input d-variate time series
  // len: Length of Time Series 
  // Short=T only computes the autocovariance up to lag  2*sqrt(len)
  
  int p = Vt.rows();
  int len = Vt.cols();
  
  Rcpp::List Gamma(len);
  Eigen::MatrixXd temp;
  for(int j=0; j<len-1; j++){
    temp = Eigen::MatrixXd::Zero(p,p);
    
    for(int k=0; k<len-j-1; k++){
      temp = temp + Vt.col(k) * Vt.col(k+j+1).transpose() / double(len);
    }
    Gamma(j) = temp;
    
  } // j=1:len-1
  
  // special case for j<-0, it is placed at Gamma(len) i.e Gamma(T)
  temp = Eigen::MatrixXd::Zero(p,p);
  for(int k=0; k<len; k++){
    temp = temp + Vt.col(k) * Vt.col(k).transpose() / double(len);
  }
  Gamma(len-1) = temp;
  
  return Gamma;
}


// [[Rcpp::export]]
Eigen::MatrixXd EvalGammaJC(Rcpp::List Gamma, int j, int len)
{
  
  // A function to evaluate Gamma(j) 
  // for any value of j from negative through zero to positive integer
  // Gamma<-Computed Gamma series for j<-0 to T-1 its dimension is (d,d,T)
  // j: index of Gamma
  // len: Length of Time Series 
  
  Eigen::MatrixXd GJ;
  if(j==0){
    GJ = Gamma(len-1);
  }else if(j<0){
    int k = -j;
    Eigen::MatrixXd temp = Gamma(k-1);
    GJ = temp.transpose();
  }else{
    GJ = Gamma(j-1);
  }
  
  return GJ;
}

// [[Rcpp::export]]
Rcpp::List CmpRhoC(Rcpp::List Gamma, int len)
{
  
  Rcpp::List Rho(len);
  Eigen::MatrixXd Gam0 = EvalGammaJC(Gamma, 0, len);
  Eigen::MatrixXd norm = Gam0.diagonal() * Gam0.diagonal().transpose();
  Eigen::MatrixXd normsqr = norm.array().pow(0.5);
  
  Eigen::MatrixXd temp;
  for(int m=0; m<len-1; m++){
    temp = Gamma(m);
    Rho(m) = temp.array() / normsqr.array();
  }
  temp = Gamma(len-1);
  Rho(len-1) = temp.array() / normsqr.array();
  
  return Rho;
}


// [[Rcpp::export]]
Eigen::MatrixXd EvalRhoMC(Rcpp::List Rho, int m, int len)
{
  // A function to evaluate Rho(m) 
  // for any value of nonnegative integer m 
  // Rho<-Computed correlogram/cross-correlogram for m<-0 to T-1 its dimension is (d,d,T)
  // m<-index of Rho
  // len<-T<-Length of Time Series 
  
  Eigen::MatrixXd RJ;
  if(m==0){
    RJ = Rho(len-1);
  }else{
    RJ = Rho(m-1);
  }
  
  return RJ;
}


// [[Rcpp::export]]
Eigen::MatrixXd CmpHatSC(Rcpp::List Rho, double C0, int KT, double cef, int p, int len)
{
  // A function to compute Bandwidth Matrix HatS from Rho
  // Rho<-Estimated correlogram/cross-correlogram Matrix of dimension (p,p,T)
  // d<-Dimension of each input vector V(t)
  // len<-T<-Length of Time Series 
  // cef<-cef value given
  
  Eigen::MatrixXd HatS = Eigen::MatrixXd::Zero(p, p);
  Eigen::MatrixXd HatQ = Eigen::MatrixXd::Ones(p, p) * 10000;
  
  int LEN = len - KT - 1;
  double val = C0*sqrt(log10(len)/len);
  int found;
  
  for(int j=0; j<p; j++){
    
    for(int k=0; k<p; k++){
      
      for(int q=1; q<LEN+1; q++){
        
        found = 1;
        
        for(int m=0; m<KT+1; m++){
          
          Eigen::MatrixXd temp = EvalRhoMC(Rho, q+m, len);
          Eigen::MatrixXd RhoJKQM = temp.array().abs();
          
          if(RhoJKQM(j,k) >= val){
            found = 0;
            break;
          }
        } // m
        
        if(found==1){
          HatQ(j,k) = double(q);
          break;
        }
      } // q
      
      if(HatQ(j,k)==10000){
        HatQ(j,k) = floor(sqrt(len));
      }
      
    }  // k
  }  // j
  
  for(int j=0; j<p; j++){
    for(int k=0; k<p; k++){
      double x;
      
      if(j==k){
        x = max(ceil(HatQ(j,k)/cef), double(1));
        HatS(j,k) = x;
      }else{
        double hatq = max(HatQ(j,k), HatQ(k,j));
        x = max(ceil(hatq/cef), double(1));
        HatS(j,k) = HatS(k,j) = x;
      }
      
    } // k
  } // j
  
  return HatS;
}


// ===================================================================
//  compute spectral density estimates
// ===================================================================
// [[Rcpp::export]]
Rcpp::List SpecEstC(Rcpp::List Gamma, int n, int p, int r, int K, Eigen::MatrixXd cross_indices, 
                    Eigen::VectorXd J_set, double l_band, double flag_c)
{
  Rcpp::List x_spc(K);
  // std::complex<double> i(0,1);
  
  // compute cross spectral
  Eigen::MatrixXcd newspec;
  for(int k=0; k<K; k++)
  {
    newspec = MatrixXcd::Zero(p,p);
    for(int h=-n; h<n+1; h++)
    {
      //cout << h << endl;
      double  q = TaperFlatC(h/l_band, flag_c);
      Eigen::MatrixXd GammaM = EvalGammaJC(Gamma, h, n);
      newspec.real() = newspec.real() + q*cos(double(h)*J_set(k))*GammaM;
      newspec.imag() = newspec.imag() - q*sin(double(h)*J_set(k))*GammaM;
      
    } // h
    x_spc(k) = newspec / (2*M_PI);
  } // k
  
  return x_spc;
}

// ===================================================================
//  compute C matrix
// ===================================================================
// [[Rcpp::export]]
Rcpp::List CEst2C(Eigen::MatrixXd x, Rcpp::List Gamma, int n_tilde, int n, int p, int r,
                  Eigen::MatrixXd cross_indices, int l_band)
{
  Rcpp::List Chat(n_tilde);
  Eigen::MatrixXd GammaM = EvalGammaJC(Gamma, n, n);
  
  for(int t=0; t<n_tilde; t++)
  {
    Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(r,2*l_band+1);
    for(int j=0; j<r; j++)
    {
      int k1 = cross_indices(j,0)-1;
      int k2 = cross_indices(j,1)-1;
      temp(j,l_band) = (x(k1,t+l_band) * x(k2,t+l_band)-GammaM(k1,k2)) / (2*M_PI);
      //if(t==0){cout << x(k1,t+l_band) * x(k2,t+l_band)-GammaM(k1,k2) << endl;}
      
      
      for(int q=0; q<l_band; q++)
      {
        Eigen::MatrixXd GammaJ1 = EvalGammaJC(Gamma,l_band-q,n);
        temp(j,q) = (x(k1,t+q) * x(k2,t+l_band) - GammaJ1(k2,k1)) / (2*M_PI);
        Eigen::MatrixXd GammaJ2 = EvalGammaJC(Gamma,q+1,n);
        temp(j,q+l_band+1) = (x(k1,t+l_band+q+1) * x(k2,t+l_band) - GammaJ2(k1,k2)) / (2*M_PI);
      } // q
    } // j
    
    Chat(t) = temp;
  } // t
  
  return Chat;
}


// ===================================================================
//  compute kernel matrix
// ===================================================================
// [[Rcpp::export]]
double BandEstC(Eigen::MatrixXd Chat, int n_tilde, int r, int l_band, int type)
{
  double a_hat,bw=0;
  Eigen::VectorXd rho = (Chat.leftCols(n_tilde-1)*Chat.rightCols(n_tilde-1).transpose()).diagonal().array()/
    (Chat.leftCols(n_tilde-1)*Chat.leftCols(n_tilde-1).transpose()).diagonal().array();
  Eigen::VectorXd sig = ((Chat.rightCols(n_tilde-1)-rho.asDiagonal()*Chat.leftCols(n_tilde-1)).array().square().rowwise().sum())/
    double(n_tilde-1);
  
  double numer_sum = 0.0;
  double denom_sum = 0.0;
  
  if(type==1){
    for(int i=0; i<(r*(2*l_band+1)); i++){
      numer_sum += 4*rho(i)*rho(i)*pow(sig(i),2)*pow((1-rho(i)),-8);
      denom_sum += pow(sig(i),2)*pow((1-rho(i)),-4);
    }
    a_hat = numer_sum/denom_sum;
    bw = 1.3221*pow((n_tilde*a_hat),0.2);
  }
  else if(type==2){
    for(int i=0; i<(r*(2*l_band+1)); i++){
      numer_sum += 4*rho(i)*rho(i)*pow(sig(i),2)*pow((1-rho(i)),-6)*pow((1+rho(i)),-2);
      denom_sum += pow(sig(i),2)*pow((1-rho(i)),-4);
    }
    a_hat = numer_sum/denom_sum;
    bw = 1.1447*pow((n_tilde*a_hat),0.33333333);
  }
  return bw;
}

// [[Rcpp::export]]
Eigen::MatrixXd etaC(int n, int p, int B, int n_tilde, double bn, int type){
  Eigen::MatrixXd kenel = Eigen::MatrixXd::Ones(n_tilde, n_tilde);
  if(type==1){
    for(int i=0; i<n_tilde; i++){
      for(int j=0; j<n_tilde; j++){
        if(i!=j){
          double temp = double(i-j)/bn;
          kenel(i,j) = 25/(12*M_PI*M_PI*temp*temp) * (sin(6*M_PI*temp/5)/(6*M_PI*temp/5) -
            cos(6*M_PI*temp/5));
        }
      }
    }
  }
  if(type==2){
    for(int i=0; i<n_tilde; i++){
      for(int j=0; j<n_tilde; j++){
        double temp = abs(double(i-j)/bn);
        if(temp<=1)
          kenel(i,j) = 1-temp;
        else
          kenel(i,j) = 0;
      }
    }
  }
  static default_random_engine e(time(0));
  static normal_distribution<double> normal(0.0,1.0);
  Eigen::MatrixXd Xi_temp=Eigen::MatrixXd::Zero(B,n_tilde);
  for(int i=0; i<B; i++){
    for(int j=0; j<n_tilde; j++){
      Xi_temp(i,j) = normal(e);    
    }
  }
  Eigen::EigenSolver<Eigen::MatrixXd> eig(kenel);
  Eigen::VectorXd EigenValue =  eig.eigenvalues().real().array();
  for(int i=0;i<n_tilde;i++){
    if(EigenValue(i)<0)EigenValue(i)=log(double(p))/double(n);
  }
  EigenValue = EigenValue.array().sqrt();
  Eigen::MatrixXd D = EigenValue.asDiagonal();
  Eigen::MatrixXd EigenVector = eig.eigenvectors().real();
  Eigen::MatrixXd Xi = EigenVector*D*EigenVector.transpose()*Xi_temp.transpose();
  return Xi;
}


// ===================================================================
//  compute long-run covariance matrix
// ===================================================================

// [[Rcpp::export]]
Eigen::MatrixXd LongCovEstC(int n_tilde, int ln, int r, Eigen::VectorXi Shat_c, 
                            Eigen::MatrixXd Chat, int Kern)
{
  Eigen::MatrixXd SigmaMat = Eigen::MatrixXd::Zero(r*(2*ln+1), r*(2*ln+1));
  for(int j1=0; j1<r*(2*ln+1); j1++)
  {
    for(int j2=0; j2<(j1+1); j2++)
    {
      int bn = Shat_c(j1,j2);
      
      if(Kern==1)
      {
        // kernel.mat is from qs kernel with b.band bandwidth
        double sig = 0.0;
        for(int q=0; q<(2*n_tilde-1); q++){
          int q_tmp = q - n_tilde + 1;
          double Kernel = TaperQsC( double(q_tmp)/bn );
          Eigen::MatrixXd pi12;
          if(q<n_tilde){
            pi12 = (Chat.block(j2,0,1,q+1) * 
              Chat.block(j1,n_tilde-q-1,1,q+1).transpose()).array() / n_tilde;
          }else{
            pi12 = (Chat.block(j1,0,1,2*n_tilde-q-1) * 
              Chat.block(j2,q-n_tilde+1,1,2*n_tilde-q-1).transpose()).array() / n_tilde;
          }
          sig += pi12(0,0)*Kernel;
          
        } // q=1:(2n-1)
        
        SigmaMat(j1, j2) = sig;
        SigmaMat(j2, j1) = sig;
        
      } // kern==1
      
      if(Kern==2){
        // kernel.mat is from qs kernel with b.band bandwidth
        double sig = 0.0;
        for(int q=0; q<(2*n_tilde-1); q++){
          int q_tmp = q - n_tilde + 1;
          double Kernel = TaperBartC( double(q_tmp)/bn );
          Eigen::MatrixXd pi12;
          if(q<n_tilde){
            pi12 = (Chat.block(j2,0,1,q+1) * 
              Chat.block(j1,n_tilde-q-1,1,q+1).transpose()).array() / n_tilde;
          }else{
            pi12 = (Chat.block(j1,0,1,2*n_tilde-q-1) * 
              Chat.block(j2,q-n_tilde+1,1,2*n_tilde-q-1).transpose()).array() / n_tilde;
          }
          sig += pi12(0,0)*Kernel;
          
        } // q=1:(2n-1)
        
        SigmaMat(j1, j2) = sig;
        SigmaMat(j2, j1) = sig;
        
      } // kern==2
      
    } // j2
  } // j1
  
  return SigmaMat;
}

