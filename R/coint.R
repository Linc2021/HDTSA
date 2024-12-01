#' @name Coint
#' @title Identifying the cointegration rank of nonstationary vector time series
#'
#' @description \code{Coint()} deals with cointegration analysis for high-dimensional
#' vector time series proposed in Zhang, Robinson and Yao (2019). Consider the model:
#' \deqn{{\bf y}_t = {\bf Ax}_t\,,}
#' where \eqn{{\bf A}} is a \eqn{p \times p} unknown and invertible constant matrix,
#'  \eqn{{\bf x}_t = ({\bf x}'_{t,1}, {\bf x}'_{t,2})'} is a latent
#'  \eqn{p \times 1} process, \eqn{{\bf x}_{t,2}} is an \eqn{r \times 1} \eqn{I(0)} process,
#'  \eqn{{\bf x}_{t,1}} is a process with nonstationary components, and no linear
#'  combination of \eqn{{\bf x}_{t,1}} is \eqn{I(0)}. This function aims to estimate the
#'  cointegration rank \eqn{r} and the invertible constant matrix \eqn{{\bf A}}.
#' 
#' @details
#' Write \eqn{\hat{\bf x}_t=\hat{\bf A}'{\bf y}_t\equiv (\hat{x}_t^1,\ldots,\hat{x}_t^p)'}.
#' When \code{type = "acf"}, \code{Coint()} estimates \eqn{r} by
#'  \deqn{\hat{r}=\sum_{i=1}^{p}1\bigg\{\frac{S_i(m)}{m}<c_0 \bigg\}} for some
#'  constant \eqn{c_0\in (0,1)} and some large constant \eqn{m}, where
#' \eqn{S_i(m)} is the sum of the sample autocorrelations of
#' \eqn{\hat{x}^{i}_{t}} over lags 1 to \eqn{m},
#' which is specified in Section 2.3 of Zhang, Robinson and Yao (2019).
#' 
#' When \code{type = "urtest"}, \code{Coint()} estimates \eqn{r} by unit root
#' tests. For \eqn{i= 1,\ldots, p}, consider the null hypothesis 
#' \deqn{H_{0,i}:\hat{x}_t^{p-i+1} \sim I(0)\,.} The estimation procedure for
#' \eqn{r} can be implemented as follows:
#' 
#' \emph{Step 1}. Start with \eqn{i=1}. Perform the unit root test proposed
#' in Chang, Cheng and Yao (2021) for \eqn{H_{0,i}}.
#' 
#' \emph{Step 2}. If the null hypothesis is not rejected at the significance
#' level \eqn{\alpha}, increment \eqn{i} by 1 and repeat Step 1. Otherwise, stop
#' the procedure and denote the value of \eqn{i} at termination as \eqn{i_0}.
#' The cointegration rank is then estimated as \eqn{\hat{r}=i_0-1}.
#' 
#' 
#' 
#'   
#' @param Y An \eqn{n \times p} data matrix \eqn{{\bf Y} = ({\bf y}_1, \dots , {\bf y}_n )'},
#'  where \eqn{n} is the number of the 
#'   observations of the \eqn{p \times 1} time series \eqn{\{{\bf y}_t\}_{t=1}^n}.
#' @param lag.k The time lag \eqn{K} used to calculate the nonnegative definte
#'   matrix \eqn{\hat{{\bf W}}_y}: \deqn{\hat{\mathbf{W}}_y\ =\
#'   \sum_{k=0}^{K}\hat{\mathbf{\Sigma}}_y(k)\hat{\mathbf{\Sigma}}_y(k)'\,,}
#'    where \eqn{\hat{\bf \Sigma}_y(k)} is the sample autocovariance of
#'   \eqn{ {\bf y}_t} at lag \eqn{k}. The default is 5.
#' @param type The method used to identify the cointegration rank. Available
#' options include: \code{"acf"} (the default) for the method based on the sample
#' autocorrelations, \code{"urtest"} for the method based on the unit root tests,
#' and \code{"both"} to apply these two methods. See Section 2.3 of Zhang, Robinson
#' and Yao (2019) and 'Details' for more information.
#' @param c0 The prescribed constant \eqn{c_0} involved in the method based on
#' the sample correlations, which is used
#' when \code{type = "acf"} or \code{type = "both"}. See Section 2.3 of Zhang, Robinson
#' and Yao (2019) and 'Details' for more information. The default is 0.3.
#' @param m The prescribed constant \eqn{m} involved in the method based on
#' the sample correlations, which is used
#' when \code{type = "acf"} or \code{type = "both"}. See Section 2.3 of Zhang, Robinson
#' and Yao (2019) and 'Details' for more information. The default is 20.
#' @param alpha The significance level \eqn{\alpha} of the unit root tests,
#' which is used when \code{type = "urtest"} or \code{type = "both"}.
#' See 'Details'. The default is 0.01.

#' @return An object of class \code{"coint"}, which contains the following
#'   components:
#'   \item{A}{The estimated \eqn{\hat{\bf A}}. }
#'   \item{coint_rank}{The estimated cointegration rank \eqn{\hat{r}}.}
#'   \item{lag.k}{The time lag used in function.}
#'   \item{method}{A string indicating which method is used to identify the
#'   cointegration rank.}
#'
#' @references 
#' Chang, J., Cheng, G., & Yao, Q. (2022).  Testing for unit
#' roots based on sample autocovariances. \emph{Biometrika}, \strong{109}, 543--550.
#' \doi{doi:10.1093/biomet/asab034}.
#' 
#' Zhang, R., Robinson, P., & Yao, Q. (2019). Identifying cointegration by
#' eigenanalysis. \emph{Journal of the American Statistical Association},
#' \strong{114}, 916--927. \doi{doi:10.1080/01621459.2018.1458620}.
#' @export
#' @importFrom stats cor
#' @useDynLib HDTSA
#' @examples
#' # Example 1 (Example 1 in Zhang, Robinson and Yao (2019))
#' ## Generate yt
#' p <- 10
#' n <- 1000
#' r <- 3
#' d <- 1
#' X <- mat.or.vec(p, n)
#' X[1,] <- arima.sim(n-d, model = list(order=c(0, d, 0)))
#' for(i in 2:3)X[i,] <- rnorm(n)
#' for(i in 4:(r+1)) X[i, ] <- arima.sim(model = list(ar = 0.5), n)
#' for(i in (r+2):p) X[i, ] <- arima.sim(n = (n-d), model = list(order=c(1, d, 1), ar=0.6, ma=0.8))
#' M1 <- matrix(c(1, 1, 0, 1/2, 0, 1, 0, 1, 0), ncol = 3, byrow = TRUE)
#' A <- matrix(runif(p*p, -3, 3), ncol = p)
#' A[1:3,1:3] <- M1
#' Y <- t(A%*%X)
#' 
#' Coint(Y, type = "both")
Coint <- function(Y, lag.k=5, type=c("acf", "urtest", "both"),
                  c0 = 0.3, m = 20, alpha = 0.01){
  
  Y <- as.matrix(Y)
  type <- match.arg(type)
  n <- nrow(Y)
  p <- ncol(Y)
  Y <- t(Y)
  
  #step1 caculate the W
  W <- diag(rep(0,p))
  mean_y <- as.matrix(rowMeans(Y))
  storage.mode(mean_y) <- "double"
  for(k in 0:lag.k){
    Sigma_y <- sigmak(Y, mean_y, k, n)
    W <- W + MatMult(Sigma_y, t(Sigma_y))
  }
  #step2 Eigendecoposition
  t <- eigen(W, symmetric=T)
  ev <- t$value
  G <- as.matrix(t$vectors)
  
  # t = svd(V0+S)
  # ev = c(t$d);
  # G = as.matrix(t$u)
  
  if (ev[p]<0)
  {cat("some eigenvalues are false", "\n")
  }
  
  ##### c2. compute the covariance and determine the cointegration rank####
  
  Z=MatMult(t(G), Y)
  Z=t(Z)
  
  ### indentify r_hat, two way:
  if (type == "acf" || type == "both"){
    z <- c(rep(0, m))
    b <- c(rep(0, p))
    
    for (i in 1:p)
    { 
      for(k in 1:m)
        z[k] <- (cor(Z[(k+1):n, i], Z[1:(n-k), i]))
      b[i] <- sum(z[1:m])
    }
    r_hat1 = sum(b< m*c0)
    if (type == "acf"){
      names(r_hat1) <- "The estimated number of cointegration rank"
      names(lag.k) <-"Time lag"
      return(structure(list(A = G, coint_rank = r_hat1, lag.k = lag.k,
                            method = "acf method"),
                       class = "coint"))
    }
  }
  # # (d) Using PP-test with p value 0.01
  # if (type == "pptest" || type == "all")
  # { 
  #   v=c(rep(0, p))
  #   for (h in p:1)
  #   {
  #     v[h] <- PP.test(Z[ , h], lshort = FALSE)$p.value
  #     if(v[h] > alpha) break
  #   }
  #   r_hat2 = p-h
  #   
  #   if (type == "pptest"){
  #     METHOD <- c("Identifying the cointegration rank of given time series",
  #                 "Using Phillips-Perron test")
  #     names(r_hat2) <- "The estimated number of cointegration rank"
  #     names(lag.k) <- "Time lag"
  #     return(structure(list(A = G, coint_rank = r_hat2, lag.k = lag.k,
  #                           method = METHOD),
  #                      class = "coint"))
  #   }
  # }
  # (d) Using unit-root test based on sample covariance with p value 0.01
  if (type == "urtest" || type == "both")
  {
    v=c(rep(0, p))
    for (l in p:1)
    {
      v[l] <- UR_test(Z[ , l], lagk.vec = 2, alpha = alpha)$reject[1, 1]
      if(v[l] == 1) {break}
    }
    r_hat2 = p - l
    if (type == "urtest"){
      names(r_hat2) <- "The estimated number of cointegration rank"
      names(lag.k) <- "Time lag"
      return(structure(list(A = G, coint_rank = r_hat2, lag.k = lag.k,
                            method = "unit root test"),
                       class = "coint"))
    }
  }
  if(type == "both"){
    result <- t(as.matrix(c(r_hat1, r_hat2)))
    rownames(result) <- "r_hat"
    colnames(result) <- c("acf", "urtest")
    names(lag.k) <-"Time lag"
    return(structure(list(A = G, coint_rank = result, lag.k = lag.k,
                          method = "both two methods"),
                     class = "coint"))
  }
}
