#' @name WN_test
#' @title Testing for white noise hypothesis in high dimension
#' @description \code{WN_test()} is the test proposed in Chang, Yao and Zhou
#' (2017) for the following hypothesis testing problems: \deqn{H_0:\{{\bf x}_t
#' \}_{t=1}^n\mathrm{\ is\ white\ noise\ \ versus\ \ }H_1:\{{\bf x}_t
#' \}_{t=1}^n\mathrm{\ is\ not\ white\ noise.} }
#'
#' @param X \eqn{{\bf X} = \{{\bf x}_1, \dots , {\bf x}_n \}'}, an \eqn{n\times
#'   p} sample matrix, where \eqn{n} is the sample size and \eqn{p} is the
#'   dimension of \eqn{{\bf x}_t}.
#' @param lag.k Time lag \eqn{K}, a positive integer, used to calculate the test
#'   statistic [See (4) in Chang, Yao and Zhou (2017)]. Default is \code{lag.k}
#'   \eqn{=2}.
#' @param B Bootstrap times for generating multivariate normal distributed
#'   random vectors in calculating the critical value. Default is \code{B}
#'   \eqn{=2000}.
#' @param pre Logical value which determines whether to performs preprocessing
#'   procedure on data matrix \code{X} or not, see Remark 1 in Chang, Yao and
#'   Zhou (2017) for more information. If \code{TRUE}, then the segment
#'   procedure will be performed to data \code{X} first. The additional
#'   parameter which control the preprocessing procedure could be set 
#'   through parameter \code{control}.
#' @param alpha The prescribed significance level. Default is 0.05.
#' @param kernel.type String, an option for choosing the symmetric kernel used
#'   in the estimation of long-run covariance matrix, for example, \code{'QS'}
#'   (Quadratic spectral kernel), \code{'Par'} (Parzen kernel) and \code{'Bart'}
#'   (Bartlett kernel), see Andrews (1991) for more information. Default option
#'   is\code{kernel.type = 'QS'}.
#' @param  control A list of control parameters. See ‘Details’ and ‘Details’ in \link{PCA_TS}.
#' @seealso \code{\link{PCA_TS}}
#' 
#' @details
#' The control argument is a list that can supply any of the following 
#'   components to \code{PCA_TS()}:
#'   \itemize{
#'   \item \code{k0}: A positive integer specified to calculate \eqn{\widehat{{\bf
#'   W}}_y}. See parameter \code{lag.k} in \code{\link{PCA_TS}} for more
#'   information.
#'   \item \code{thresh}: Logical. It determines whether to perform the threshold method
#'   to estimate \eqn{\widehat{{\bf W}}_y} or not. See parameter \code{thresh}
#'   in \code{\link{PCA_TS}} for more information.
#'   \item \code{tuning.vec}: The value of thresholding tuning parameter \eqn{\lambda}.
#'   See parameter \code{tuning.vec} in \code{\link{PCA_TS}} for more
#'   information.
#'   \item \code{standardize}: Whether the variables will be standardized to
#'   have mean zero and unit standard deviation. Default \code{FALSE}.
#'   \item \code{opt}: Method options for calculating ransformation matrix in preprocessing
#'   procedure on data matrix \code{X}. See parameter details in \link{PCA_TS}.
#'   Default options is \code{opt=1}.
#'   \item \code{K}: The number of folders used in the cross validation for the
#'   selection of \eqn{\lambda}, the default is 5. It is required when
#'   \code{control(thresh = TRUE)}. See parameter details in \link{PCA_TS}.
#'   }
#' 
#'
#' @return An object of class "hdtstest" is a list containing the following
#'   components:
#'
#'   \item{statistic}{The value of the test statistic.}
#'   \item{p.value}{Numerical value which represents the p-value of the test
#'   based on the observed data \eqn{\{{\bf x}_t\}_{t=1}^n}.}
#'   \item{lag.k}{The time lag used in function.}
#'   \item{method}{A character string indicating what method was performed.}
#'   \item{kernel.type}{A character string indicating what kenel method was performed.}
#'   
#' @references Chang, J., Yao, Q., & Zhou, W. (2017). \emph{Testing for
#'   high-dimensional white noise using maximum cross-correlations}, Biometrika,
#'   Vol. 104, pp. 111–127.
#'
#'   Chang, J., Guo, B., & Yao, Q. (2018). \emph{Principal component analysis for
#'   second-order stationary vector time series}, The Annals of Statistics, Vol.
#'   46, pp. 2094–2124.
#'
#'   Cai, T., and Liu, W. (2011). \emph{Adaptive thresholding for sparse
#'   covariance matrix estimation},  Journal of the American Statistical
#'   Association, Vol. 106, pp. 672--684.
#'
#' @examples
#' n <- 200
#' p <- 10
#' X <- matrix(rnorm(n*p),n,p)
#' res <- WN_test(X)
#' Pvalue <- res$p.value
#' rej <- res$reject
#' @useDynLib HDTSA
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @export
# WN_test = function(X, lag.k = 2, B = 1000, method = c("CYZ","CLL"),
#                    kernel.type = c('QS','Par','Bart'), resampling = FALSE,
#                    pre = FALSE, alpha = 0.05,k0 = 5, thresh = FALSE,
#                    tuning.vec = NULL,
#                    opt = 1, control = list()){
WN_test = function(X, lag.k = 2, B = 1000, kernel.type = c("QS", "Par", "Bart"),
                   pre = FALSE, alpha = 0.05, control = list()){
  
  kernel.type <- match.arg(kernel.type)
  ken_type <- switch(kernel.type,
                     "QS" = 1,
                     "Par" = 2,
                     "Bart" = 3)
  
  if (pre == TRUE){
    con <- list(k0 = 5, thresh = FALSE, tuning.vec = NULL, K=5, opt = 1)
    con[(namc <- names(control))] <- control
    # print(con)
    X_pre <- segmentTS(X, lag.k = con$k0,
                      thresh = con$thresh,
                      tuning.vec = con$tuning.vec,
                      opt = con$opt,
                      control = control,
                      K = con$K)
    
    X <- X_pre$Z
    Bx <- X_pre$B
  }
  n <- nrow(X)
  p <- ncol(X)
  
  Tn_list <- WN_teststatC(X,n,p,lag.k)
  Tn <- Tn_list$Tn
  sigma_zero <- Tn_list$sigma_zero
  X_mean <- Tn_list$X_mean
  
  ft <- WN_ftC(n, lag.k, p, X, X_mean)
  bn <- bandwith(ft, lag.k, p, p, ken_type)
  
  Gnstar <- WN_bootc(n, lag.k, p, B, bn, ken_type, ft, X, sigma_zero) # critical value
  p.value <- mean(Gnstar > Tn)
  # Results = list(reject = (p.value<0.05), p.value = p.value)
  
  names(Tn) <- "Statistic"
  names(lag.k) <-"Time lag"
  names(kernel.type) <- "Symmetric kernel"
  METHOD <- "Testing for white noise hypothesis in high dimension"
  structure(list(statistic = Tn, p.value = p.value, lag.k=lag.k,
                 method = METHOD, 
                 kernel = kernel.type),
            class = "hdtstest")
}
