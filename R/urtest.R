#' Testing for unit roots based on sample autocovariances
#'
#' The test proposed in Chang, Cheng and Yao (2021) for the following hypothesis
#' testing problems: \deqn{H_0:Y_t \sim I(0)\ \ \mathrm{versus}\ \ H_1:Y_t \sim
#' I(d)\ \mathrm{for\ some\ integer\ }d \geq 2.}
#' @param Y \eqn{Y = \{y_1, \dots , y_n \}}, the observations of a univariate
#'   time series used for the test.
#' @param lagk.vec Time lag \eqn{K_0} used to calculate the test statistic, see
#'   Section 2.1 in Chang, Cheng and Yao (2021). It can be a vector containing
#'   more than one time lag. If it is a vector, the procedure will output all
#'   the test results based on the different \eqn{K_0} in the vector
#'   \code{lagk.vec}. If \code{lagk.vec} is missing, the default value we choose
#'   lagk.vec=c(0,1,2,3,4).
#' @param con_vec Constant \eqn{c_\kappa}, see (5) in Chang, Cheng and Yao
#'   (2021). It also can be a vector. If missing, the default value we use 0.55.
#' @param alpha The prescribed significance level. Default is 0.05.

#' @return A dataframe containing the following components:
#'
#'   \item{result}{\code{'1'} means we reject the null hypothesis and \code{'0'}
#'   means we do not reject the null hypothesis.}
#'
#' @references Chang, J., Cheng, G. & Yao, Q. (2021).  \emph{Testing for unit
#'   roots based on sample autocovariances}. Available at
#'   \url{https://arxiv.org/abs/2006.07551}
#' @export
#' @importFrom sandwich lrvar
#' @importFrom stats lm
#' @useDynLib HDTSA
#' @importFrom Rcpp sourceCpp
#' @importFrom stats qnorm
#' @importFrom Rcpp evalCpp
#' @import Rcpp
#' @examples
#' N=100
#' Y=arima.sim(list(ar=c(0.9)), n = 2*N, sd=sqrt(1))
#' con_vec=c(0.45,0.55,0.65)
#' lagk.vec=c(0,1,2)
#' ur.test(Y,lagk.vec=lagk.vec, con_vec=con_vec,alpha=0.05)
#' ur.test(Y,alpha=0.05)


ur.test <- function(Y, lagk.vec=lagk.vec, con_vec=con_vec, alpha=alpha) {
  
  args = as.list(match.call())
  if(is.null(args$lagk.vec)){
    lagk.vec=c(0,1,2,3,4);
  }
  args = as.list(match.call())
  if(is.null(args$con_vec)){
    con_vec=0.55;
  }
  Tnvec=NULL; nm=NULL; colnm=NULL;
  
  for (i in con_vec)colnm = c(colnm, paste("con=", i, sep=""))
  
  for(kk in 1:length(lagk.vec)){
    
    K0=lagk.vec[kk]+1                       #eg. K0=1, gamma(0)
    nm=c(nm,paste("K0=", K0-1, sep=""))
    
    
    n=length(Y)                           ## sample size
    N=floor(n/2)
    N1=2*N-K0
    
    sgn_matrix=matrix(0,N1,K0)            ### sign matrix
    for(t in 1:N1){
      for(k in 1:K0){
        sgn_matrix[t,k]=sign(k+t-N-1-0.5) ##  eg. K0=1, gamma(0)
      }
    }
    
    Y=Y;  DY=diff(Y) ## diffential Y
    
    au_Y =drop(acf(Y,lag.max =K0+1, type = c("covariance"),plot = FALSE)$acf)   ## gamma(Y)
    au_DY=drop(acf(DY,lag.max=K0+1, type = c("covariance"),plot = FALSE)$acf)   ## gamma(X)
    
    short_Var=var(DY)                                         ## shortrun variance
    long_Var=n*lrvar(DY, type = "andrews", prewhite = FALSE)  ## longrun variance
    ratio_Var=short_Var/long_Var                              ## variance ratio
    
    
    ## esttimate rho
    Z2=Y[1:(n-1)];  Z1=Y[2:n]
    DZ2=diff(Z2);   DZ1=diff(Z1)
    rho_hat=lm(DZ2~DZ1)$coefficients[2]; bb=1+rho_hat         ## rho_hat
    
    au_Ratio=(au_Y[1]+au_Y[2])/(au_DY[1]+au_DY[2])            ## ratio
    
    Y1=Y[1:N]                                   ##data spliting
    Y2=Y[(N+1):(2*N)]
    
    ## auto covariance
    auto_cov=drop(acf(Y,lag.max =K0+1,type = c("covariance"),plot = FALSE)$acf)
    auto_cov1=drop(acf(Y1,lag.max =K0+1,type = c("covariance"),plot = FALSE)$acf)
    auto_cov2=drop(acf(Y2,lag.max =K0+1,type = c("covariance"),plot = FALSE)$acf)
    
    T1=sum((auto_cov1[1:(K0)])^2)
    T2=sum((auto_cov2[1:(K0)])^2)             ## test statistics
    
    ## construct Qt
    N1=2*N-K0; Y=Y-mean(Y)
    ft=Y[1:N1]%*%t(rep(1,K0));
    for(t in 1:N1){
      ft[t,]=ft[t,]*Y[(t):(t+K0-1)]                       ## gamma(0), .... gamma(K0-1) data
    }
    gamma_hat=t(auto_cov[1:(K0)]%*%t(rep(1,N1)))         ## estimate gamma(0),.... gamma(K0-1)
    ytk=2*(ft-gamma_hat)*sgn_matrix                       ## construct ytk
    sgn_Auto=t(sign(auto_cov1[1:(K0)])%*%t(rep(1,N1)))    ## sign auto covariance
    xitk=2*ytk*gamma_hat; Qt=apply(xitk,1,sum)               ## Qt
    lr_Qt=lrvar(Qt, type = "andrews", prewhite = FALSE)   ## long-run variance Qt
    
    ## test procedure
    
    kappa=2/(ratio_Var*bb);                          ## kappa
    if(lr_Qt>0){                                     ## variance >0
      #cv=qnorm(1-alpha)*sqrt(lr_Qt)+T1
      #Tnvec=c(Tnvec,T2>cv)                         ## no truncated
      for(tt in 1:length(con_vec)){
        ck=con_vec[tt]
        if (au_Ratio<=(ck*kappa*N^{3/5})){
          th_d=10^5                                  ## truncated belongs to H0
        }  else{th_d=0.1*log(N)}                     ## truncated belongs to H1
        
        cv=min(qnorm(1-alpha)*sqrt(lr_Qt)+T1, th_d)
        Tnvec=c(Tnvec, T2>cv)
      }
    }
  }
  res.table=matrix(as.numeric(Tnvec), length(lagk.vec), byrow=T)
  rownames(res.table)=nm        #rownames ("K0=1", "K0=2")
  colnames(res.table) = colnm
  return(list(result=res.table))
  
}
