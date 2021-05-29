#' @name PCA4_TS
#' @title Principal component analysis for time serise
#' @description
#' \code{PCA4_TS()} seek for a contemporaneous linear transformation for
#' a multivariate time series such that the transformed series is segmented
#' into several lower-dimensional subseries, and those subseries are
#' uncorrelated with each other both contemporaneously and serially, then permutation step can
#' grouping the components of a multivariate series X into \eqn{q} groups.
#'
#'
#'
#' @details
#'  Before segment step, the time series need to be normalized then the variance of \eqn{Y}_{t} will be \eqn{I}_{p}. 
#'  When \eqn{p} is large (\eqn{>n}),it is necessary to use the package \code{clime}. When \eqn{p} is large (\eqn{p=O(n^{1/2}})), it is necessary to use the thresholding method, see more information in Chang et al. (2018).
#'
#'
#'
#'
#' @param Y  A data matrix with \eqn{n} rows and \eqn{p} columns, where \eqn{n} is the sample size and \eqn{p} is the dimension of the time series.
#' @param lag.k  A positive integer specified to calculate Wy. See (2.5) in Chang et al. (2014).
#' @param thresh   Logical. If \code{FALSE} (the default), no thresholding will be applied. If \code{TRUE}, a
#'                 thresholding method will be applied first to estimate Wy, see (3.4) and (3.5) in Chang et al. (2018).
#' @param tuning.vec  The value of thresholding parameter \eqn{\lambda}. The thresholding level is specified by
#'                  \deqn{ u = \lambda {(log p/n)}^{1/2}.}
#'                    Default value is 2. If \code{tuning.vec} is a vector, then a cross validation method proposed in Cai and Liu (2011) will be used
#'                    to choose the best tuning parameter.
#'
#' @param K   The number of folders used in the cross validation, the default is 5. It is required when \code{thresh} is \code{TRUE}.
#' @param isvol Logical. If \code{FALSE} (the default), then prewhiten each series by fitting a univariate AR model with
#'          the order between 0 and 5 determined by AIC. If \code{TRUE}, then prewhiten each volatility process using GARCH(1,1) model.
#' @param permutation The method of permutation procdure,the option is \code{'max'} or \code{'fdr'},default option is \code{'max'}
#' @param m A positive constant used to calculate the maximum cross correlation over the lags between \eqn{-m} and \eqn{m}. If \eqn{m} is not specified, the default constant \eqn{10*log10(n/p)}
#'          will be used.
#' @param Beta The error rate in FDR.
#' @param just4pre Logical. If \code{TRUE}, then only segment procdure will be used, If \code{FALSE} (the default), the nomal procdure will be used.
#' @export
#' @return The first step of segment containing the following components:
#' \item{B}{The \eqn{p*p} transformation matrix such that \eqn{x_t = By_t}}
#' \item{X}{The transformed series with \eqn{n} rows and \eqn{p} columns}
#' @return  Output of the second step is a list containing the following components:
#' \item{NoGroups}{Number of groups with at least two components series}
#' \item{Nos_of_Members}{Number of members in each of groups with at least two members}
#' \item{Groups}{Indices of components in each of groups with at least two members}
#' \item{maxcorr}{Maximum correlation (over lags) of \eqn{p(p-1)/2} pairs in descending order if \code{permutation} is \code{max}}
#' \item{corrRatio}{Ratios of successive values from maxcorr if \code{permutation} is \code{max}}
#' \item{Pvalues}{Pvalue for multiple test H_0 for each of \eqn{p(p-1)/2} pairs in ascending order if \code{permutation} is \code{fdr}}
#' \item{NoConnectedPairs}{Number of connected pairs}
#' \item{Xpre}{The prewhitened data with \eqn{n-R} rows and \eqn{p} columns}
#'
#' @seealso \code{\link{segmentTS}}, \code{\link{permutationMax}}, \code{\link{permutationFDR}}.
#' @note The first step is transform the time series. It calculate linear transformation of the \eqn{p}-variate time series (or volatility processes) \eqn{y_t} such that the transformed series \eqn{x_t=By_t} is segmented into several
#' lower-dimensional subseries, and those subseries are uncorrelated with each other both contemporaneously and serially.
#' The second step is grouping the transformed time series, The permutation is determined by grouping the components of a multivariate series X into \eqn{q} groups, where \eqn{q} and the cardinal numbers of those groups are also unknown. See also \code{\link{permutationFDR}} and \code{\link{permutationMax}}
#'
#'
#' @references Chang, J., Guo, B. & Yao, Q. (2018). \emph{Principal component analysis for second-order stationary vector time series}. The Annals of Statistics, Vol. 46, pp. 2094–2124.
#'
#'             Cai, T. and Liu, W. (2011). \emph{Adaptive thresholding for sparse covariance matrix estimation}.  Journal of the American Statistical Association 106: 672-684.
#'
#'             Cai, T.T., Liu, W., and Luo, X. (2011). \emph{A constrained `1 minimization approach for sparse precision matrix estimation}. Journal of the American Statistical Association 106(494): 594-607.
## @family aggregate functions
#'

## Perform Step 2 (i.e. permutation) of Chang, Guo and Yao (2014)
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
#

PCA4_TS <- function(Y, lag.k=5, thresh=FALSE, tuning.vec=NULL, K=5, 
                    isvol=FALSE, permutation=c('max',"fdr"), m=NULL, Beta=NULL, 
                    just4pre=FALSE)
{
  #for timeseries
  if(is.logical(just4pre)){
    print("Wrong parameter type (just4pre)")
    print("Set default option just4pre=FALSE")
  }
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
	  output=list(NoGroups=out$NoGroups, Nos_of_Members=out$Nos_of_Members, Groups=out$Groups, maxcorr=out$maxcorr,
	              corrRatio=out$corrRatio,NoConnectedPairs=out$r,Xpre=out$Xpre,B=seglist$B)
	return(output)
  }
  	#permutation of FDR
	else if(permutation=='fdr')
	  {
      out=permutationFDR(X,isvol, Beta, m)
      output=list(NoGroups=out$NoGroups, Nos_of_Members=out$Nos_of_Members, Groups=out$Groups, Pvalues=out$Pvalues,
                  NoConnectedPairs=out$r,Xpre=out$Xpre,B=seglist$B)
		return(output)
      }
}


