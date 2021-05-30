#' @name WN_test
#' @title Testing for white noise hypothesis in high dimension
#' @description
#' \code{WN_test()} is a new omnibus test for vector white noise using the maximum 
#' absolute auto-correlations and cross-correlations of the component series. \deqn{H_0:\{{\bf x}_t \}_{t=1}^n\mathrm{\ is\ white\ noise\ \ versus\ \ }H_1:\{{\bf x}_t \}_{t=1}^n\mathrm{\ is\ not\ white\ noise.} }
#' Based on an approximation by the \eqn{L_\infty} norm of a normal random vector, the critical 
#' value of the test can be evaluated by bootstrapping from a multivariate normal distribution.
#' @param X \eqn{{\bf X} = \{{\bf x}_1, \dots , {\bf x}_n \}'}, a \eqn{n\times p} data matrix used testing whether the data are white noise process.
#' @param lag.k Time lag \eqn{K}, a positive integer, used to calculate test statistic. Default is \code{lag.k} \eqn{=2}.
#' @param B Bootstrap times for generate vector from a multivariate normal distribution. Default is \code{B} \eqn{=2000}
#' @param pre Logical which determines whether performs preprocessing procedure on data matrix \eqn{{\bf X}}. The segment procedure will performs to data when parameter \code{just4pre} in \code{\link{PCA4_TS}} set to true. The 
#' three additional options including \code{thresh}, \code{tuning.vec} and \code{cv.num} are the same as 
#' those in \code{\link{PCA4_TS}}. This makes the maximum correlation greater, and therefore the test is more powerful when \eqn{p} is large or moderately large.
#' @param alpha Significance level used for testing. Default is 0.05.
#' @param kernel.type An option for choosing an optimal symmetric kernel type, for example \code{'QS'}, \code{'Par'}, \code{'Bart'}. Default is \code{kernel.type='QS'}.          
#' @param k0 A positive integer specified to calculate \eqn{\widehat{{\bf W}}_y}. See parameter \code{lag.k} in \code{\link{PCA4_TS}} for more information.
#' @param thresh Logical. It determines whether to perform the threshold method to estimate \eqn{\widehat{{\bf W}}_y}. See parameter \code{thresh} in \code{\link{PCA4_TS}} for more information.
#' @param tuning.vec The value of thresholding parameter \eqn{\lambda}. See parameter \code{tuning.vec} in \code{\link{PCA4_TS}} for more information.
#' @seealso \code{\link{PCA4_TS}}
#' 
#' @return An object of class "WN_test" is a list containing the following components:
#'
#' \item{reject}{Logical value which represents whether the data are white noise process.}
#' \item{p.value}{Numerical value which represents the \eqn{p}-value obtained in the bootstrap sampling process}
#' @references Chang, J., Yao, Q. & Zhou, W. (2017). \emph{Testing for high-dimensional white noise using maximum cross-correlations}. Biometrika, Vol. 104, pp. 111–127.
#' 
#'             Chang, J., Guo, B. & Yao, Q. (2018). \emph{Principal component analysis for second-order stationary vector time series}. The Annals of Statistics, Vol. 46, pp. 2094–2124.
#'
#'             Cai, T. and Liu, W. (2011). \emph{Adaptive thresholding for sparse covariance matrix estimation}.  Journal of the American Statistical Association 106: 672-684.
#'            
#' @examples 
#' n <- 200
#' p <- 150
#' X <- matrix(rnorm(n*p),n,p)
#' res <- WN_test(X)
#' Pvalue <- res$p.value
#' rej <- res$reject
#' @useDynLib HDTSA
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @import Rcpp 
#' @import RcppEigen
#' @export

WN_test = function(X, lag.k=2, B=2000, kernel.type='QS', pre=FALSE,alpha=0.05,k0=5,thresh=FALSE,
                   tuning.vec=NULL){
  
  if (pre==T){
    X_pre = PCA4_TS(X, lag.k=k0, thresh=thresh, just4pre = pre, tuning.vec = tuning.vec)
    #X_pre = SgmentTS(X, k0); #PCA4_TS(X, k0, thresh, just4pre = TRUE, tuning.vec = tuningvec)
    X = X_pre$X
    Bx = X_pre$B
  }
  n = nrow(X)
  p = ncol(X)
  
  Tn_list = WN_teststatC(X,n,p,lag.k)
  Tn = Tn_list$Tn
  sigma_zero = Tn_list$sigma_zero
  X_mean = Tn_list$X_mean
  
  if(kernel.type =="QS")ken_type=1
  if(kernel.type =="PR")ken_type=2
  if(kernel.type =="BT")ken_type=3
  ft = WN_ftC(n,lag.k,p,X,X_mean)
  bn = WN_bandwith(ft,lag.k,p,ken_type)
  
  Gnstar = WN_bootc(n,lag.k,p,B,bn,ken_type,ft,X,sigma_zero) # critical value
  p.value = mean(Gnstar>Tn)
  Results = list(reject = (p.value<0.05), p.value = p.value)
  
  return(Results)
}
