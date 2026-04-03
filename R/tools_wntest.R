#' @name WN_test
#' @title Testing for white noise hypothesis in high dimension
#' @description \code{WN_test()} implements the white noise tests proposed in Chang, Yao and Zhou
#' (2017) and Chang et al. (2026+) for the following hypothesis testing problem: \deqn{H_0:\{{\bf y}_t
#' \}_{t=1}^n\mathrm{\ is\ white\ noise\ \ versus\ \ }H_1:\{{\bf y}_t
#' \}_{t=1}^n\mathrm{\ is\ not\ white\ noise.} }
#'
#' @param Y An \eqn{n \times p} data matrix \eqn{{\bf Y} = ({\bf y}_1, \dots , {\bf y}_n )'},
#' where \eqn{n} is the number of the observations of the \eqn{p \times 1}
#' time series \eqn{\{{\bf y}_t\}_{t=1}^n}.
#' @param lag.k The time lag \eqn{K} used to calculate the test
#'   statistic [See (4) of Chang, Yao and Zhou (2017) and (6) of 
#'   Chang et al. (2026+).]. The default is 2.
#' @param B The number of bootstrap replications for calculating the critical value.
#' The default is 1000.
#' @param pre Logical. This parameter is used when \code{method = "L_inf"}.
#' If \code{TRUE} (the default), the time series PCA
#' proposed in Chang, Guo and Yao (2018) should be performed on
#' \eqn{\{{\bf y}_t\}_{t=1}^n} before implementing the white noise test [See Remark 1
#' of Chang, Yao and Zhou (2017)]. The time series PCA is implemented by using
#' the function \code{\link{PCA_TS}} with the arguments passed by \code{control.PCA}.
#' 
#' @param alpha The significance level of the test. The default is 0.05.
#' @param method The option for method used in the white noise test. Available options
#' include: \code{"L_inf"} (the default) for the method proposed by Chang, Yao and
#' Zhou (2017), and \code{"L_2"} for the method proposed by Chang et al. (2026+).

#' @param kernel.type The option for choosing the symmetric kernel used
#'   in the estimation of long-run covariance matrix, which is used when
#'   \code{method = "L_inf"}. Available options include:
#'   \code{"QS"} (the default) for the Quadratic spectral kernel, \code{"Par"}
#'   for the Parzen kernel, and \code{"Bart"} for the Bartlett kernel.
#'   See Chang, Yao and Zhou (2017) for more information. 
#'   
#' @param control.PCA A list of control arguments passed to the function
#' \code{PCA_TS} when \code{method = "L_inf"} and \code{pre = TRUE}, including
#' \code{lag.k}, \code{opt}, \code{thresh}, \code{delta}, and the associated
#' arguments passed to the \code{clime} function
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
#'   Chang, J., He, J., Li, W., & Lin, C. (2026+). An adaptive \eqn{L_2}-type test
#'   for high-dimensional white noise. \emph{Preprint}.
#'   
#'   Chang, J., Guo, B., & Yao, Q. (2018). Principal component analysis for
#'   second-order stationary vector time series. \emph{The Annals of Statistics},
#'   \strong{46}, 2094--2124. \doi{doi:10.1214/17-AOS1613}.
#'   
#'   Chang, J., Yao, Q., & Zhou, W. (2017). Testing for
#'   high-dimensional white noise using maximum cross-correlations. \emph{Biometrika},
#'   \strong{104}, 111--127. \doi{doi:10.1093/biomet/asw066}.
#'   
#'
#' @examples
#' #Example 1
#' ## Generate data
#' n <- 200
#' p <- 10
#' Y <- matrix(rnorm(n * p), n, p)
#' 
#' ## L_inf
#' res1 <- WN_test(Y, method ="L_inf")
#' Pvalue1 <- res1$p.value
#' statistic1 <- res1$statistic
#' 
#' ## L_2
#' res2 <- WN_test(Y, method = "L_2")
#' Pvalue2 <- res2$p.value
#' statistic2 <- res2$statistic
#' 
#' @useDynLib HDTSA
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @export
#' 
# todo method = c("L_inf","L_2")
WN_test = function(Y, lag.k = 2, B = 1000, method = c("L_inf","L_2"),
                   kernel.type = c("QS", "Par", "Bart"),
                   pre = FALSE, alpha = 0.05, control.PCA = list()){
  method <- match.arg(method)
  n <- nrow(Y)
  p <- ncol(Y)
  if(method == 'L_inf'){
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
  else if(method == 'L_2'){
    n <- nrow(Y)
    p <- ncol(Y)
    Hn <- rep(0, lag.k)
    const <- sqrt(2) / (n - 1)
    tmp <- Y %*% t(Y)
    diag(tmp) <- 0
    sigma_n1 <- const * sum(tmp^2)
    for(lag in c(1:lag.k)){
      Gn_1 <- sum(tmp * tmp[c(c((n-lag+1):n),c(1:(n-lag))),c(c((n-lag+1):n),c(1:(n-lag)))])
      if(lag > 1){
        Hn[lag] <- Hn[lag-1] + Gn_1/sigma_n1
      }
      else{
        Hn[lag] <- Gn_1/sigma_n1
      }
    }
    #--------------------calculate cv for Hn----------------------------------
    
    Hn_B <- resampling(Y,n,p,B,lag.k)
    cv_vec_resampling <- Hn_B[,floor(B*(1-alpha))]
    
    names(Hn) <- sprintf("Statistic (time lag = %d)", 1:lag.k)
    names(lag.k) <-"Time lag"
    p.value = rowMeans(sweep(Hn_B, 1, Hn, FUN = ">"))
    
    return(structure(list(statistic = Hn, p.value = p.value, lag.k=lag.k),
              class = "hdtstest"))

  }
}
