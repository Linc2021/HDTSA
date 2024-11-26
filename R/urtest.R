#' @name UR_test
#' @title Testing for unit roots based on sample autocovariances
#'
#' @description The test implements the test proposed in Chang, Cheng and Yao (2022)
#' for the following hypothesis testing problems:
#' \deqn{H_0:Y_t \sim I(0)\ \ \mathrm{versus}\ \ H_1:Y_t \sim
#' I(d)\ \mathrm{for\ some\ integer\ }d \geq 2\,,} where \eqn{\{Y_t\}_{t=1}^n} is
#' a univariate time series.
#' @param Y A vector \eqn{Y = (y_1, \dots , y_n )}, where \eqn{n} is the number
#' of the observations.
#' @param lagk.vec A vector of length \eqn{d} with each element representing a
#' specified time lag \eqn{K_0} used to calculate the test statistic
#' [See Section 2.1 of Chang, Cheng and Yao (2022)].
#' The procedure outputs all test results corresponding to the specified \eqn{d} lags.
#' The default is \code{c(0, 1, 2, 3, 4)}.
#' @param con_vec A vector of length \eqn{m} with each element representing a
#' specified constant \eqn{c_\kappa} [See (5) of Chang, Cheng and Yao (2022)]. 
#' The procedure outputs all test results corresponding to \eqn{m} values of \eqn{c_\kappa}.
#' The default is \code{c(0.55)}.
#' @param alpha The significance level of the test. Default is 0.05.

#' @return An object of class \code{"urtest"}, which contains the following
#'   components:
#'
#'   \item{statistic}{A vector of length \eqn{d} with each element representing
#'   the value of a test statistic corresponding to one of the \eqn{d} time lags.}
#'   \item{reject}{An \eqn{m \times d} data matrix containing test results for
#'   various parameter (\eqn{c_\kappa,K_0}) combinations. Each row
#'    corresponds to a specific value of \eqn{c_\kappa} and 
#'    each column represents a specific time lag \eqn{K_0},  A value of
#'    \code{'1'} indicates rejection of the null hypothesis,
#'    while \code{'0'} indicates non-rejection.}
#'   \item{lag.vec}{The time lag used in function.}
#'   
#'
#' @references Chang, J., Cheng, G., & Yao, Q. (2022). Testing for Unit
#'   Roots Based on Sample Autocovariances. \emph{Biometrika}, \strong{109},
#'   543--550. \doi{doi:10.1093/biomet/asab034}
#'   
#' @export
#' @importFrom sandwich lrvar
#' @importFrom stats lm
#' @useDynLib HDTSA
#' @importFrom Rcpp sourceCpp
#' @importFrom stats qnorm
#' @importFrom Rcpp evalCpp
#' @examples
#' # Example 1
#' ## generate yt
#' N <- 100
#' Y <-arima.sim(list(ar = c(0.9)), n = 2*N, sd = sqrt(1))
#' con_vec <- c(0.45, 0.55, 0.65)
#' lagk.vec <- c(0, 1, 2)
#' 
#' UR_test(Y, lagk.vec = lagk.vec, con_vec = con_vec, alpha = 0.05)
#' UR_test(Y, alpha = 0.05)


UR_test <- function(Y, lagk.vec = NULL, con_vec = NULL, alpha = 0.05) {
  
  args <- as.list(match.call())
  if(is.null(args$lagk.vec)){
    lagk.vec <- c(0,1,2,3,4);
  }
  args <- as.list(match.call())
  if(is.null(args$con_vec)){
    con_vec <- 0.55;
  }
  Tnvec <- NULL; nm <- NULL; colnm <- NULL; statistic_vec <- NULL
  
  for (i in con_vec)colnm <- c(colnm, paste("con=", i, sep = ""))
  
  for(kk in 1:length(lagk.vec)){
    
    K0 <- lagk.vec[kk]+1                       #eg. K0=1, gamma(0)
    nm <- c(nm, paste("time_lag=", K0-1, sep = ""))
    
    
    n <- length(Y)                           ## sample size
    N <- floor(n/2)
    N1 <- 2*N - K0
    
    sgn_matrix <- matrix(0, N1, K0)            ### sign matrix
    for(t in 1:N1){
      for(k in 1:K0){
        sgn_matrix[t, k] <- sign(k+t-N-1-0.5) ##  eg. K0=1, gamma(0)
      }
    }
    
    Y <- Y;  DY <- diff(Y) ## diffential Y
    
    au_Y <- drop(acf(Y, lag.max = K0+1, type = c("covariance"), plot = FALSE)$acf)   ## gamma(Y)
    au_DY <- drop(acf(DY, lag.max = K0+1, type = c("covariance"), plot = FALSE)$acf)   ## gamma(X)
    
    short_Var <- var(DY)                                         ## shortrun variance
    long_Var <- n * lrvar(DY, type = "andrews", prewhite = FALSE)  ## longrun variance
    ratio_Var <- short_Var / long_Var                              ## variance ratio
    
    
    ## esttimate rho
    Z2 <- Y[1:(n-1)];  Z1 <- Y[2:n]
    DZ2 <- diff(Z2);   DZ1 <- diff(Z1)
    rho_hat <- lm(DZ2~DZ1)$coefficients[2]; bb <- 1 + rho_hat         ## rho_hat
    
    au_Ratio <- (au_Y[1] + au_Y[2]) / (au_DY[1] + au_DY[2])            ## ratio
    
    Y1 <- Y[1:N]                                   ##data spliting
    Y2 <- Y[(N+1):(2*N)]
    
    ## auto covariance
    auto_cov <- drop(acf(Y, lag.max = K0+1, type = c("covariance"), plot = FALSE)$acf)
    auto_cov1 <- drop(acf(Y1, lag.max = K0+1, type = c("covariance"), plot = FALSE)$acf)
    auto_cov2 <- drop(acf(Y2, lag.max = K0+1, type = c("covariance"), plot = FALSE)$acf)
    
    T1 <- sum((auto_cov1[1:(K0)])^2)
    T2 <- sum((auto_cov2[1:(K0)])^2)             ## test statistics
    
    ## construct Qt
    N1 <- 2 * N - K0; Y <- Y - mean(Y)
    ft <- Y[1:N1] %*% t(rep(1, K0));
    for(t in 1:N1){
      ft[t, ] <- ft[t, ] * Y[(t):(t+K0-1)]                       ## gamma(0), .... gamma(K0-1) data
    }
    gamma_hat <- t(auto_cov[1:(K0)] %*% t(rep(1, N1)))         ## estimate gamma(0),.... gamma(K0-1)
    ytk <- 2 * (ft - gamma_hat) * sgn_matrix                       ## construct ytk
    sgn_Auto <- t(sign(auto_cov1[1:(K0)]) %*% t(rep(1, N1)))    ## sign auto covariance
    xitk <- 2 * ytk * gamma_hat; Qt <- apply(xitk, 1, sum)               ## Qt
    lr_Qt <- lrvar(Qt, type = "andrews", prewhite = FALSE)   ## long-run variance Qt
    
    ## test procedure
    
    kappa <- 2 / (ratio_Var*bb);                          ## kappa
    if(lr_Qt > 0){                                     ## variance >0
      #cv=qnorm(1-alpha)*sqrt(lr_Qt)+T1
      #Tnvec=c(Tnvec,T2>cv)                         ## no truncated
      for(tt in 1:length(con_vec)){
        ck <- con_vec[tt]
        if (au_Ratio <= (ck * kappa * N^{3/5})){
          th_d <- 10^5                                  ## truncated belongs to H0
        }  else{th_d <- 0.1 * log(N)}                     ## truncated belongs to H1
        
        cv <- min(qnorm(1-alpha) * sqrt(lr_Qt) + T1, th_d)
        statistic_vec <- c(statistic_vec, T2)
        Tnvec <- c(Tnvec, T2>cv)
      }
    }
  }
  res.table <- matrix(as.numeric(Tnvec), length(lagk.vec), byrow=T)
  rownames(res.table) <- nm        #rownames ("K0=1", "K0=2")
  colnames(res.table) <- colnm
  # statistic_vec <- matrix(as.numeric(statistic_vec), length(lagk.vec))
  # colnames(statistic_vec) <- "statistic"
  # METHOD <- "Testing for unit roots based on sample autocovariances"
  # return(list(result=res.table))
  structure(list(statistic = statistic_vec, reject = res.table,
                 lag.k = lagk.vec), class = "urtest")
  
}
