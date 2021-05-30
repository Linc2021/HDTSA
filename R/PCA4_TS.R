#' @name PCA4_TS
#' @title Principal component analysis for time serise
#' @description
#' \code{PCA4_TS()} seek for a contemporaneous linear transformation for
#' a multivariate time series such that the transformed series is segmented
#' into several lower-dimensional subseries: \deqn{{\bf y}_t={\bf Ax}_t,} where \eqn{{\bf x}_t} is an unobservable \eqn{p \times 1} weakly stationary time series consisting of \eqn{q\ (\geq 1)} both contemporaneously and serially
#'  uncorrelated subseries. then permutation step can grouping the components of a multivariate series \eqn{{\bf x}_t} into \eqn{q} groups.
#'
#'
#'
#' @details
#'  Before segment step, the time series \eqn{{\bf y}_t} need to be normalized then the variance of \eqn{{\bf y}_t} will be \eqn{{\bf I}_p}. 
#'  When \eqn{p} is large (\eqn{>n}), it is necessary to use the package \pkg{clime} to estimate the precision matrix \eqn{\widehat{{\bf V}}^{-1}} where \eqn{\widehat{{\bf V}}} is a consistent estimator for \eqn{\mathrm{Var}({\bf y}_t)}. In segment step, when \eqn{p} is large (\eqn{p=O(n^{1/2}})), it is recommended to use the thresholding method to calculate \eqn{\widehat{{\bf W}}_y}, see more information in Chang et al. (2018).
#'
#'
#'
#'
#' @param Y  \eqn{{\bf Y} = \{{\bf y}_1, \dots , {\bf y}_n \}'}, a data matrix with \eqn{n} rows and \eqn{p} columns, where \eqn{n} is the sample size and \eqn{p} is the dimension of the time series.
#' @param lag.k  Time lag \eqn{k_0} used to calculate the nonnegative definte matrix \eqn{\widehat{{\bf W}}_y}: \deqn{\widehat{\mathbf{W}}_y\ =\ \sum_{k=0}^{k_0}\widehat{\mathbf{\Sigma}}_y(k)\widehat{\mathbf{\Sigma}}_y(k)'=\mathbf{I}_p+\sum_{k=1}^{k_0}\widehat{\mathbf{\Sigma}}_y(k)\widehat{\mathbf{\Sigma}}_y(k)', }
#'               where \eqn{\widehat{\bf \Sigma}_y(k)} is the sample autocovariance of \eqn{ {\bf y}_t} at lag \eqn{k}. See (2.5) in Chang et al. (2018).
#' @param thresh   Logical. If \code{FALSE} (the default), no thresholding will be applied. If \code{TRUE}, a
#'                 thresholding method will be applied first to estimate \eqn{\widehat{{\bf W}}_y}, see (3.3) and (3.5) in Chang et al. (2018).
#' @param tuning.vec  The value of thresholding parameter \eqn{\lambda}. The thresholding level is specified by
#'                  \deqn{ u = \lambda {(log p/n)}^{1/2},}
#'                    where default value is 2. If \code{tuning.vec} is a vector, then a cross validation method proposed in Cai and Liu (2011) will be used
#'                    to choose the best tuning parameter \eqn{\lambda}.
#'
#' @param K   The number of folders used in the cross validation, the default is 5. It is required when \code{thresh} is \code{TRUE}.
#' @param isvol Logical. If \code{FALSE} (the default), then prewhiten each series by fitting a univariate AR model with
#'          the order between 0 and 5 determined by AIC. If \code{TRUE}, then prewhiten each volatility process using GARCH(1,1) model.
#' @param permutation The method of permutation procdure, where option is \code{'max'} or \code{'fdr'}, default is \code{permutation='max'}
#' @param m A positive constant used to calculate the maximum cross correlation over the lags between \eqn{-m} and \eqn{m}. If \eqn{m} is not specified, then default constant \eqn{10*log10(n/p)}
#'          will be used.
#' @param Beta The error rate in FDR.
#' @param just4pre Logical. If \code{TRUE}, then only segment procdure will be used, if \code{FALSE} (the default), then the nomal procdure (including permutation step) will be used. See \code{\link{WN_test}} for more application
#' @export
#' @return The first step output of segment containing the following components:
#' \item{B}{The \eqn{p\times p} transformation matrix such that \eqn{{\bf x}_t = {\bf By}_t}}
#' \item{X}{The transformed series with \eqn{n} rows and \eqn{p} columns}
#' @return  The second step output is a list containing the following components:
#' \item{NoGroups}{Number of groups with at least two components series}
#' \item{Nos_of_Members}{Number of members in each of groups with at least two members}
#' \item{Groups}{Indices of components in each of groups with at least two members}
#' \item{maxcorr}{Maximum correlation (over lags) of \eqn{p(p-1)/2} pairs in descending order if \code{permutation = 'max'}}
#' \item{corrRatio}{Ratios of successive values from maxcorr if \code{permutation= 'max'}}
#' \item{Pvalues}{Pvalue for multiple test: \deqn{H_0:\rho_{i,j}(h)=0\ \mathrm{\ for\ any}\ h\ =\ 0,\pm 1 ,\pm2,\dots,\pm m} for each of \eqn{p(p-1)/2} pairs in ascending order if \code{permutation = 'fdr'}, see more information in Chang et al. (2018).}
#' \item{NoConnectPairs}{Number of connected pairs}
#' \item{Xpre}{The prewhitened data with \eqn{n-R} rows and \eqn{p} columns}
#'
#' @note The first step is transform the time series. It calculate linear transformation of the \eqn{p}-variate time series (or volatility processes) \eqn{y_t} such that the transformed series \eqn{x_t=By_t} is segmented into several
#' lower-dimensional subseries.
#' The second step is grouping the transformed time series, The permutation is determined by grouping the components of a multivariate series X into \eqn{q} groups, where \eqn{q} and the cardinal numbers of those groups are also unknown.
#'
#'
#' @references Chang, J., Guo, B. & Yao, Q. (2018). \emph{Principal component analysis for second-order stationary vector time series}. The Annals of Statistics, Vol. 46, pp. 2094–2124.
#'
#'             Cai, T. and Liu, W. (2011). \emph{Adaptive thresholding for sparse covariance matrix estimation}.  Journal of the American Statistical Association 106: 672-684.
#'
#'             Cai, T.T., Liu, W., and Luo, X. (2011). \emph{A constrained `1 minimization approach for sparse precision matrix estimation}. Journal of the American Statistical Association 106(494): 594-607.
## @family aggregate functions
#'

## Perform Step 2 (i.e. permutation) of Chang, Guo and Yao (2018)
## using the maximum cross correlation method
##
## Output1: $NoGroups -- No. of groups with at least two components series
## Output2: $NosOfMembersInGroups -- number of members in each of groups
## 				     with at least two members
## Output3: $Groups -- indices of components in each of groups with at least two members
## Output4: $Pvalues -- Pvalue for multiple test H_0 for each of p(p-1)/2 pairs
##          in ascending order
## Output5: $NoConnectedPairs -- No of connected pairs by FDR at rate Beta
#' @importFrom stats acf ar pnorm var 
#' @useDynLib HDTSA
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @import Rcpp 
#' @import RcppEigen
#' @export
#' @examples 
#' ## Example 1 (Example 5 of Chang et al.(2018)). 
#' ## p=6, x_t consists of 3 independent subseries with 3, 2 and 1 components.    
#'
#' p=6;n=1500
#' # Generate x_t
#' X=mat.or.vec(p,n)
#' x=arima.sim(model=list(ar=c(0.5, 0.3), ma=c(-0.9, 0.3, 1.2,1.3)),n=n+2,sd=1)
#' for(i in 1:3) X[i,]=x[i:(n+i-1)]  
#' x=arima.sim(model=list(ar=c(0.8,-0.5),ma=c(1,0.8,1.8) ),n=n+1,sd=1)
#' for(i in 4:5) X[i,]=x[(i-3):(n+i-4)]   
#' x=arima.sim(model=list(ar=c(-0.7, -0.5), ma=c(-1, -0.8)),n=n,sd=1)
#' X[6,]=x
#' # Generate y_t 
#' A=matrix(runif(p*p, -3, 3), ncol=p)
#' Y=A%*%X  
#' Y=t(Y)
#' res=PCA4_TS(Y, lag.k=5,permutation = "max")
#' # The transformed series z_t 
#' Z=res$X 
#' # Plot the cross correlogram of z_t and y_t
#' Y=data.frame(Y);Z=data.frame(Z)
#' names(Y)=c("Y1","Y2","Y3","Y4","Y5","Y6")
#' names(Z)=c("Z1","Z2","Z3","Z4","Z5","Z6")
#' # The cross correlogram of y_t shows no block pattern 
#' acfY=acf(Y) 
#' # The cross correlogram of z_t shows 3-2-1 block pattern  
#' acfZ=acf(Z)
#'       
#' ## Example 2 (Example 6 of Chang et al.(2018)).
#' ## p=20, x_t consists of 5 independent subseries with 6, 5, 4, 3 and 2 components.    
#'
#' p=20;n=3000
#' # Generate x_t
#' X=mat.or.vec(p,n)
#' x=arima.sim(model=list(ar=c(0.5, 0.3), ma=c(-0.9, 0.3, 1.2,1.3)),n.start=500,n=n+5,sd=1)
#' for(i in 1:6) X[i,]=x[i:(n+i-1)]
#' x=arima.sim(model=list(ar=c(-0.4,0.5),ma=c(1,0.8,1.5,1.8)),n.start=500,n=n+4,sd=1)
#' for(i in 7:11) X[i,]=x[(i-6):(n+i-7)]
#' x=arima.sim(model=list(ar=c(0.85,-0.3),ma=c(1,0.5,1.2)), n.start=500,n=n+3,sd=1)
#' for(i in 12:15) X[i,]=x[(i-11):(n+i-12)]
#' x=arima.sim(model=list(ar=c(0.8,-0.5),ma=c(1,0.8,1.8)),n.start=500,n=n+2,sd=1)
#' for(i in 16:18) X[i,]=x[(i-15):(n+i-16)]
#' x=arima.sim(model=list(ar=c(-0.7, -0.5), ma=c(-1, -0.8)),n.start=500,n=n+1,sd=1)
#' for(i in 19:20) X[i,]=x[(i-18):(n+i-19)]
#' # Generate y_t 
#' A=matrix(runif(p*p, -3, 3), ncol=p)
#' Y=A%*%X  
#' Y=t(Y)
#' res=PCA4_TS(Y, lag.k=5,permutation = "max")
#' # The transformed series z_t 
#' Z=res$X 
#' # Plot the cross correlogram of x_t and y_t
#' Y=data.frame(Y);Z=data.frame(Z)
#' namesY=NULL;namesZ=NULL
#' for(i in 1:p)
#' {
#'    namesY=c(namesY,paste0("Y",i))
#'    namesZ=c(namesZ,paste0("Z",i))
#' }  
#' names(Y)=namesY;names(Z)=namesZ
#' # The cross correlogram of y_t shows no block pattern 
#' acfY=acf(Y, plot=FALSE)
#' plot(acfY, max.mfrow=6, xlab='', ylab='',  mar=c(1.8,1.3,1.6,0.5), 
#'      oma=c(1,1.2,1.2,1), mgp=c(0.8,0.4,0),cex.main=1)    
#' # The cross correlogram of z_t shows 6-5-4-3-2 block pattern  
#' acfZ=acf(Z, plot=FALSE)
#' plot(acfZ, max.mfrow=6, xlab='', ylab='',  mar=c(1.8,1.3,1.6,0.5),
#'      oma=c(1,1.2,1.2,1), mgp=c(0.8,0.4,0),cex.main=1)
#' # Identify the permutation mechanism
#' permutation=res
#' permutation$Groups  

PCA4_TS <- function(Y, lag.k=5, thresh=FALSE, tuning.vec=NULL, K=5, 
                    isvol=FALSE, permutation=c('max',"fdr"), m=NULL, Beta, 
                    just4pre=FALSE)
{
  #for timeseries
  permutation <- match.arg(permutation)
  
  seglist=segmentTS(Y=Y,lag.k = lag.k,isvol = isvol, 
                      thresh = thresh,
                      tuning.vec = tuning.vec,
                      K = K)
  X=seglist$X
  B=seglist$B
  if(just4pre==TRUE){
    Yt=list(B=B, X=X)
    return(Yt)
  }
  #permutation of MAX
  if(permutation == 'max')
  {
  	out=permutationMax(X, isvol, m)
	  output=list(B=B, X=X,NoGroups=out$NoGroups, Nos_of_Members=out$Nos_of_Members, Groups=out$Groups, maxcorr=out$maxcorr,
	              corrRatio=out$corrRatio,NoConnectedPairs=out$r,Xpre=out$Xpre,B=seglist$B)
	return(output)
  }
  	#permutation of FDR
	else if(permutation=='fdr')
	  {
      out=permutationFDR(X,isvol, Beta, m)
      output=list(B=B, X=X,NoGroups=out$NoGroups, Nos_of_Members=out$Nos_of_Members, Groups=out$Groups, Pvalues=out$Pvalues,
                  NoConnectPairs=out$r,Xpre=out$Xpre,B=seglist$B)
		return(output)
      }
}


