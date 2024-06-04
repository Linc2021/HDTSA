#' @name PCA_TS
#' @title Principal component analysis for time serise
#' @description \code{PCA_TS()} seeks for a contemporaneous linear
#'   transformation for a multivariate time series such that the transformed
#'   series is segmented into several lower-dimensional subseries: \deqn{{\bf
#'   y}_t={\bf Ax}_t,} where \eqn{{\bf x}_t} is an unobservable \eqn{p \times 1}
#'   weakly stationary time series consisting of \eqn{q\ (\geq 1)} both
#'   contemporaneously and serially uncorrelated subseries. See Chang, Guo and
#'   Yao (2018).
#'
#'
#'
#' @details When \eqn{p>n^{1/2}}, the procedure use package \pkg{clime} to
#'   estimate the precision matrix \eqn{\widehat{{\bf V}}^{-1}}, otherwise uses
#'   function \code{cov()} to estimate \eqn{\widehat{{\bf V}}} and calculate its
#'   inverse. When \eqn{p>n^{1/2}}, we recommend to use the thresholding method
#'   to calculate \eqn{\widehat{{\bf W}}_y}, see more information in Chang, Guo
#'   and Yao (2018).
#'
#'
#'
#'
#' @param Y  \eqn{{\bf Y} = \{{\bf y}_1, \dots , {\bf y}_n \}'}, a data matrix
#'   with \eqn{n} rows and \eqn{p} columns, where \eqn{n} is the sample size and
#'   \eqn{p} is the dimension of \eqn{{\bf y}_t}. The procedure will first
#'   normalize \eqn{{\bf y}_t} as \eqn{\widehat{{\bf V}}^{-1/2}{\bf y}_t}, where
#'   \eqn{\widehat{{\bf V}}} is an estimator for covariance of \eqn{{\bf y}_t}.
#'   See details below for the selection of \eqn{\widehat{{\bf V}}^{-1}}.
#' @param lag.k  Time lag \eqn{k_0} used to calculate the nonnegative definte
#'   matrix \eqn{\widehat{{\bf W}}_y}: \deqn{\widehat{\mathbf{W}}_y\ =\
#'   \sum_{k=0}^{k_0}\widehat{\mathbf{\Sigma}}_y(k)\widehat{\mathbf{\Sigma}}_y(k)'=\mathbf{I}_p+\sum_{k=1}^{k_0}\widehat{\mathbf{\Sigma}}_y(k)\widehat{\mathbf{\Sigma}}_y(k)',
#'    } where \eqn{\widehat{\bf \Sigma}_y(k)} is the sample autocovariance of
#'   \eqn{ \widehat{{\bf V}}^{-1/2}{\bf y}_t} at lag \eqn{k}. See (2.5) in
#'   Chang, Guo and Yao (2018).
#' @param thresh   Logical. If \code{FALSE} (the default), no thresholding will
#'   be applied to estimate \eqn{\widehat{{\bf W}}_y}. If \code{TRUE}, a
#'   thresholding method will be applied first to estimate \eqn{\widehat{{\bf
#'   W}}_y}, see (3.5) in Chang, Guo and Yao (2018).
#' @param tuning.vec  The value of the tuning parameter \eqn{\lambda} in the
#'   thresholding level \eqn{ u = \lambda \sqrt{n^{-1}\log p}}, where default
#'   value is 2. If \code{tuning.vec} is a vector, then a cross validation
#'   method proposed in Cai and Liu (2011) will be used to choose the best
#'   tuning parameter \eqn{\lambda}.
#'
#' @param K   The number of folders used in the cross validation for the
#'   selection of \eqn{\lambda}, the default is 5. It is required when
#'   \code{thresh = TRUE}.
#' @param prewhiten Logical. If \code{TRUE} (the default), we prewhiten each
#'   transformed component series of \eqn{\hat{\bf z}_t} [See Section 2.2.1 in
#'   Chang, Guo and Yao (2018)] by fitting a univariate AR model with the order
#'   between 0 and 5 determined by AIC. If \code{FALSE}, then prewhiten
#'   procedure will not be performed to \eqn{\hat{\bf z}_t}.
#' @param permutation The method of permutation procedure to assign the
#'   components of \eqn{\hat{\bf z}_t} to different groups [See Section 2.2.1 in
#'   Chang, Guo and Yao (2018)]. Option is \code{'max'} (Maximum cross
#'   correlation method) or \code{'fdr'} (False discovery rate procedure based
#'   on multiple tests), default is \code{permutation = 'max'}. See Sections
#'   2.2.2 and 2.2.3 in Chang, Guo and Yao (2018) for more information.
#' @param m A positive constant used in the permutation procedure [See (2.10) in
#'   Chang, Guo and Yao (2018)]. If \eqn{m} is not specified, then default
#'   option is \code{m = }10.
#' @param beta The error rate used in the permutation procedure when
#'   \code{permutation = 'fdr'}.
#' @param just4pre Logical. If \code{TRUE}, the procedure outputs \eqn{\hat{\bf
#'   z}_t}, otherwise outputs \eqn{\hat{\bf x}_t} (the permutated version of
#'   \eqn{\hat{\bf z}_t}).
#' @param verbose Logical. If \code{TRUE}, the main results of the permutation 
#'   procedure will be output on the console. Otherwise, the result will not be 
#'   output.
#' @export
#' @return The output of the segment procedure is a list containing the
#'   following components: 
#'   \item{B}{The \eqn{p\times p} transformation matrix
#'   such that \eqn{\hat{\bf z}_t = \widehat{\bf B}{\bf y}_t}, where
#'   \eqn{\widehat{\bf B}=\widehat{\bf \Gamma}_y\widehat{{\bf V}}^{-1/2}}.}
#'   \item{Z}{\eqn{\hat{\bf Z}=\{\hat{\bf z}_1,\dots,\hat{\bf z}_n\}'}, the
#'   transformed series with \eqn{n} rows and \eqn{p} columns.}
#'   
#' @return  The output of the permutation procedure is a list containing the
#'   following components: 
#'   \item{NoGroups}{number of groups with at least two components series.}
#'   \item{No_of_Members}{The cardinalities of different groups.}
#'   \item{Groups}{The indices of the components in \eqn{\hat{\bf z}_t} that
#'   belongs to a group.}
#'   \item{method}{a character string indicating what method was performed.}
#'   
#'
#'
#' @references Chang, J., Guo, B. & Yao, Q. (2018). \emph{Principal component
#'   analysis for second-order stationary vector time series}, The Annals of
#'   Statistics, Vol. 46, pp. 2094--2124.
#'
#'   Cai, T. & Liu, W. (2011). \emph{Adaptive thresholding for sparse covariance
#'   matrix estimation},  Journal of the American Statistical Association, Vol.
#'   106, pp. 672--684.
#'
#'   Cai, T., Liu, W., & Luo, X. (2011). \emph{A constrained l1 minimization
#'   approach for sparse precision matrix estimation}, Journal of the American
#'   Statistical Association, Vol. 106, pp. 594--607.
#' @importFrom stats acf ar pnorm var
#' @useDynLib HDTSA
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @import Rcpp
#' @export
#' @examples
#' ## Example 1 (Example 5 of Chang Guo and Yao (2018)).
#' ## p=6, x_t consists of 3 independent subseries with 3, 2 and 1 components.
#'
#' p <- 6;n <- 1500
#' # Generate x_t
#' X <- mat.or.vec(p,n)
#' x <- arima.sim(model=list(ar=c(0.5, 0.3), ma=c(-0.9, 0.3, 1.2,1.3)),
#' n=n+2,sd=1)
#' for(i in 1:3) X[i,] <- x[i:(n+i-1)]
#' x <- arima.sim(model=list(ar=c(0.8,-0.5),ma=c(1,0.8,1.8) ),n=n+1,sd=1)
#' for(i in 4:5) X[i,] <- x[(i-3):(n+i-4)]
#' x <- arima.sim(model=list(ar=c(-0.7, -0.5), ma=c(-1, -0.8)),n=n,sd=1)
#' X[6,] <- x
#' # Generate y_t
#' A <- matrix(runif(p*p, -3, 3), ncol=p)
#' Y <- A%*%X
#' Y <- t(Y)
#' res <- PCA_TS(Y, lag.k=5,permutation = "max")
#' res1=PCA_TS(Y, lag.k=5,permutation = "fdr", beta=10^(-10))
#' # The transformed series z_t
#' Z <- res$Z
#' # Plot the cross correlogram of z_t and y_t
#' Y <- data.frame(Y);Z=data.frame(Z)
#' names(Y) <- c("Y1","Y2","Y3","Y4","Y5","Y6")
#' names(Z) <- c("Z1","Z2","Z3","Z4","Z5","Z6")
#' # The cross correlogram of y_t shows no block pattern
#' acfY <- acf(Y)
#' # The cross correlogram of z_t shows 3-2-1 block pattern
#' acfZ <- acf(Z)
#'
#' ## Example 2 (Example 6 of Chang Guo and Yao (2018)).
#' ## p=20, x_t consists of 5 independent subseries with 6, 5, 4, 3 and 2 components.
#' p <- 20;n <- 3000
#' # Generate x_t
#' X <- mat.or.vec(p,n)
#' x <- arima.sim(model=list(ar=c(0.5, 0.3), ma=c(-0.9, 0.3, 1.2,1.3)),n.start=500,
#' n=n+5,sd=1)
#' for(i in 1:6) X[i,] <- x[i:(n+i-1)]
#' x <- arima.sim(model=list(ar=c(-0.4,0.5),ma=c(1,0.8,1.5,1.8)),n.start=500,n=n+4,sd=1)
#' for(i in 7:11) X[i,] <- x[(i-6):(n+i-7)]
#' x <- arima.sim(model=list(ar=c(0.85,-0.3),ma=c(1,0.5,1.2)), n.start=500,n=n+3,sd=1)
#' for(i in 12:15) X[i,] <- x[(i-11):(n+i-12)]
#' x <- arima.sim(model=list(ar=c(0.8,-0.5),ma=c(1,0.8,1.8)),n.start=500,n=n+2,sd=1)
#' for(i in 16:18) X[i,] <- x[(i-15):(n+i-16)]
#' x <- arima.sim(model=list(ar=c(-0.7, -0.5), ma=c(-1, -0.8)),n.start=500,n=n+1,sd=1)
#' for(i in 19:20) X[i,] <- x[(i-18):(n+i-19)]
#' # Generate y_t
#' A <- matrix(runif(p*p, -3, 3), ncol=p)
#' Y <- A%*%X
#' Y <- t(Y)
#' res <- PCA_TS(Y, lag.k=5,permutation = "max")
#' res1 <- PCA_TS(Y, lag.k=5,permutation = "fdr",beta=10^(-200))
#' # The transformed series z_t
#' Z <- res$Z
#' # Plot the cross correlogram of x_t and y_t
#' Y <- data.frame(Y);Z <- data.frame(Z)
#' namesY=NULL;namesZ=NULL
#' for(i in 1:p)
#' {
#'    namesY <- c(namesY,paste0("Y",i))
#'    namesZ <- c(namesZ,paste0("Z",i))
#' }
#' names(Y) <- namesY;names(Z) <- namesZ
#' # The cross correlogram of y_t shows no block pattern
#' acfY <- acf(Y, plot=FALSE)
#' plot(acfY, max.mfrow=6, xlab='', ylab='',  mar=c(1.8,1.3,1.6,0.5),
#'      oma=c(1,1.2,1.2,1), mgp=c(0.8,0.4,0),cex.main=1)
#' # The cross correlogram of z_t shows 6-5-4-3-2 block pattern
#' acfZ <- acf(Z, plot=FALSE)
#' plot(acfZ, max.mfrow=6, xlab='', ylab='',  mar=c(1.8,1.3,1.6,0.5),
#'      oma=c(1,1.2,1.2,1), mgp=c(0.8,0.4,0),cex.main=1)
#' # Identify the permutation mechanism
#' permutation <- res
#' permutation$Groups  

PCA_TS <- function(Y, lag.k=5, thresh=FALSE, tuning.vec=NULL, K=5, 
                    prewhiten=TRUE, permutation=c('max',"fdr"), m=NULL, beta, 
                    just4pre=FALSE,verbose = FALSE)
{
  #for timeseries
  permutation <- match.arg(permutation)
  
  seglist=segmentTS(Y=Y,lag.k = lag.k, 
                      thresh = thresh,
                      tuning.vec = tuning.vec,
                      K = K)
  Z=seglist$Z
  B=seglist$B
  METHOD = "Principal component analysis for time serise"
  if(just4pre==TRUE){
    METHOD = c(METHOD, "Only segment procedure")
    Yt=structure(list(B=B, Z=Z, method = METHOD),
                 class = "tspca")
    return(Yt)
  }
  #permutation of MAX
  if(permutation == 'max')
  {
    METHOD = c(METHOD, "Maximum cross correlation method")
  	out=permutationMax(Z, prewhiten, m, verbose)
	  output=structure(list(B=B, Z=Z, NoGroups=out$NoGroups,
	                        No_of_Members=out$No_of_Members, Groups=out$Groups,
	                        method = METHOD),
	                  class = "tspca")
	return(output)
  }
  	#permutation of FDR
	else if(permutation=='fdr'){
	  METHOD = c(METHOD, "FDR based on multiple tests")
	  out=permutationFDR(Z,prewhiten, beta, m, verbose)
    output=structure(list(B=B, Z=Z, NoGroups=out$NoGroups, 
                          No_of_Members=out$No_of_Members, Groups=out$Groups,
                          method = METHOD),
                     class = "tspca")
		return(output)
      }
}


