#' @name WN_test
#' @title Testing for white noise hypothesis in high dimension
#' @description \code{WN_test()} implements the test proposed in Chang, Yao and Zhou
#' (2017) for the following hypothesis testing problem: \deqn{H_0:\{{\bf y}_t
#' \}_{t=1}^n\mathrm{\ is\ white\ noise\ \ versus\ \ }H_1:\{{\bf y}_t
#' \}_{t=1}^n\mathrm{\ is\ not\ white\ noise.} }
#'
#' @param Y An \eqn{n \times p} data matrix \eqn{{\bf Y} = ({\bf y}_1, \dots , {\bf y}_n )'},
#' where \eqn{n} is the number of the observations of the \eqn{p \times 1}
#' time series \eqn{\{{\bf y}_t\}_{t=1}^n}.
#' @param lag.k The time lag \eqn{K} used to calculate the test
#'   statistic [See (4) of Chang, Yao and Zhou (2017)]. The default is 2.
#' @param B The number of bootstrap replications for generating multivariate
#' normally distributed random vectors when calculating the critical value.
#' The default is 1000.
#' @param pre Logical. If \code{TRUE} (the default), the time series PCA
#' proposed in Chang, Guo and Yao (2018) should be performed on
#' \eqn{\{{\bf y}_t\}_{t=1}^n} before implementing the white noise test [See Remark 1
#' of Chang, Yao and Zhou (2017)]. The time series PCA is implemented by using
#' the function \code{\link{PCA_TS}} with the arguments passed by \code{control.PCA}.
#' @param alpha The significance level of the test. The default is 0.05.
#' @param kernel.type The option for choosing the symmetric kernel used
#'   in the estimation of long-run covariance matrix. Available options include:
#'   \code{"QS"} (the default) for the Quadratic spectral kernel, \code{"Par"}
#'   for the Parzen kernel, and \code{"Bart"} for the Bartlett kernel.
#'   See Chang, Yao and Zhou (2017) for more information.
#'   
#' @param control.PCA A list of control arguments passed to the function
#' \code{PCA_TS()}, including \code{lag.k}, \code{opt}, \code{thresh},
#' \code{delta}, and the associated arguments passed to the \code{clime} function
#' (when \code{opt = 2}). See 'Details’ in \code{\link{PCA_TS}}.
#' 
#' @seealso \code{\link{PCA_TS}}
#' 
#' 
#'
#' @return An object of class \code{"hdtstest"}, which contains the following
#'   components:
#'
#'   \item{statistic}{The test statistic of the test.}
#'   \item{p.value}{The p-value of the test.}
#'   \item{lag.k}{The time lag used in function.}
#'   \item{kernel.type}{The kernel used in function.}
#'   
#' @references 
#'   Chang, J., Guo, B., & Yao, Q. (2018). Principal component analysis for
#'   second-order stationary vector time series. \emph{The Annals of Statistics},
#'   \strong{46}, 2094--2124. \doi{doi:10.1214/17-AOS1613}.
#'
#'   Chang, J., Yao, Q., & Zhou, W. (2017). Testing for
#'   high-dimensional white noise using maximum cross-correlations. \emph{Biometrika},
#'   \strong{104}, 111--127. \doi{doi:10.1093/biomet/asw066}.
#'
#' @examples
#' #Example 1
#' ## Generate xt
#' n <- 200
#' p <- 10
#' Y <- matrix(rnorm(n * p), n, p)
#' 
#' res <- WN_test(Y)
#' Pvalue <- res$p.value
#' rej <- res$reject
#' @useDynLib HDTSA
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @export
#' 
# todo method = c("CYZ","CLL")
WN_test = function(Y, lag.k = 2, B = 1000, kernel.type = c("QS", "Par", "Bart"),
                   pre = FALSE, alpha = 0.05, control.PCA = list()){
  
  kernel.type <- match.arg(kernel.type)
  ken_type <- switch(kernel.type,
                     "QS" = 1,
                     "Par" = 2,
                     "Bart" = 3)
  
  if (pre == TRUE){
    con <- list(lag.k = 5, thresh = FALSE, delta = 2 * sqrt(log(ncol(Y)) / nrow(Y)),
                opt = 1)
    con[(namc <- names(control.PCA))] <- control.PCA
    # print(con)
    X_pre <- segmentTS(Y, lag.k = con$lag.k,
                      thresh = con$thresh,
                      delta = con$delta,
                      opt = con$opt,
                      control = control.PCA)
    
    Y <- X_pre$Z
    By <- X_pre$B
  }
  n <- nrow(Y)
  p <- ncol(Y)
  
  Tn_list <- WN_teststatC(Y,n,p,lag.k)
  Tn <- Tn_list$Tn
  sigma_zero <- Tn_list$sigma_zero
  X_mean <- Tn_list$X_mean
  
  ft <- WN_ftC(n, lag.k, p, Y, X_mean)
  bn <- bandwith(ft, lag.k, p, p, ken_type)
  
  boot_nomal <- matrix(rnorm(B*(n-lag.k)), B, n-lag.k)
  Gnstar <- WN_bootc(n, lag.k, p, B, bn, ken_type, ft, Y, sigma_zero, boot_nomal) # critical value
  p.value <- mean(Gnstar > Tn)
  # Results = list(reject = (p.value<0.05), p.value = p.value)
  
  names(Tn) <- "Statistic"
  names(lag.k) <-"Time lag"
  names(kernel.type) <- "Symmetric kernel"

  structure(list(statistic = Tn, p.value = p.value, lag.k=lag.k,
                 kernel = kernel.type),
            class = "hdtstest")
}
