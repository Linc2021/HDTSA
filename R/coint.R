#' @name Coint
#' @title Identifying cointegration rank of given time series
#'
#' @description \code{Coint} seeks for a contemporaneous linear
#'   transformation for a multivariate time series such that we can identifying 
#'   cointegration rank from the transformed series.
#'   
#' @param Y \eqn{{\bf Y} = \{{\bf y}_1, \dots , {\bf y}_n \}'}, a data matrix
#'   with \eqn{n} rows and \eqn{p} columns, where \eqn{n} is the sample size and
#'   \eqn{p} is the dimension of \eqn{{\bf y}_t}.
#' @param lag.k Time lag \eqn{k_0} used to calculate the nonnegative definte
#'   matrix \eqn{\widehat{{\bf W}}_y}: \deqn{\widehat{\mathbf{W}}_y\ =\
#'   \sum_{k=0}^{k_0}\widehat{\mathbf{\Sigma}}_y(k)\widehat{\mathbf{\Sigma}}_y(k)'}
#'    where \eqn{\widehat{\bf \Sigma}_y(k)} is the sample autocovariance of
#'   \eqn{ \widehat{{\bf y}_t}} at lag \eqn{k}.
#' @param type The method of identifying cointegration rank after segment 
#'   procedure. Option is \code{'acf'}, \code{'all'}, \code{'chang'} or \code{'pptest'}
#'   , the latter two methods use the unit-root test method to identify the 
#'   cointegration rank, and the option \code{type = 'all'} means use all three
#'   methods to identify the cointegration rank. Default is \code{type = 'acf'}.
#'    See Sections 2.3 in Zhang, Robinson and Yao (2019) for more information.
#' @param c0 The prescribed constant for identifying 
#'   cointegration rank using \code{"acf"} method. Default is 0.3.[See (2.3) in
#'   Zhang, Robinson and Yao (2019)].
#' @param m The prescribed constant for identifying 
#'   cointegration rank using \code{"acf"} method. Default is 20. [See (2.3) in
#'   Zhang, Robinson and Yao (2019)].
#' @param alpha The prescribed significance level for identifying 
#'   cointegration rank using \code{"pptest","chang"} method. Default is 0.01.
#'   [See (2.3) in Zhang, Robinson and Yao (2019)].

#' @return An object of class "coint" is a list containing the following
#'   components:
#'   \item{Z}{The transformed series with \eqn{n} rows and \eqn{p} columns. }
#'   \item{coint_rank}{A \eqn{1 \times 1} matrix representing the cointegration rank.
#'   If \code{type = 'all'}, then return a \eqn{1 \times 3} matrix representing
#'    the cointegration rank of all three methods.}
#'   \item{lag.k}{a prescribed positive integer which means the time lags used
#'    to calculate the statistic.}
#'   \item{method}{a character string indicating which method was performed.}
#'
#' @references Zhang, R., Robinson, P. & Yao, Q. (2019).  \emph{Identifying 
#' Cointegration by Eigenanalysis}.  Journal of the American Statistical 
#' Association, Vol. 114, pp. 916--927
#' @export
#' @importFrom stats PP.test cor
#' @useDynLib HDTSA
#' @examples
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
#' Coint(Y, type = "all")
Coint <- function(Y, lag.k=5, type=c("acf","pptest","Chang","all"),
                  c0 = 0.3, m = 20, alpha = 0.01){
  # Y the observed time series with length n and dimension p
  # k0 the lag order  of the covariance in recoverying  cointegration space
  
  type <- match.arg(type)
  n <- nrow(Y)
  p <- ncol(Y)
  Y <- t(Y)
  
  # S0=cov(t(Y)); V0=S0%*%S0
  # S=matrix(0, nrow=p, ncol=p)
  # for(m in 1:lag.k)
  # {Ym=Y[,(m+1):n]; Zm=Y[,1:(n-m)]
  # Sm=cov(t(Ym), t(Zm));
  # Vm=Sm%*%t(Sm);
  # S=S+Vm
  # }
  #print(V0+S)
  #step1 caculate the W
  W <- diag(rep(0,p))
  mean_y <- as.matrix(rowMeans(Y))
  storage.mode(mean_y) <- "double"
  for(k in 0:lag.k){
    Sigma_y <- sigmak(Y, mean_y, k, n)
    W <- W + MatMult(Sigma_y, t(Sigma_y))
  }
  #print(W)
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
  if (type == "acf" || type == "all"){
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
      # result <- as.matrix(r_hat1)
      METHOD <- c("Identifying cointegration rank of given time series",
                  "using acf method")
      # rownames(result) <- "r_hat"
      # colnames(result) <- "acf"
      names(r_hat1) <- "The estimated number of cointegration rank"
      names(lag.k) <-"Time lag"
      return(structure(list(Z = Z, coint_rank = r_hat1, lag.k = lag.k,
                            method = METHOD),
                       class = "coint"))
    }
  }
  # (d) Using PP-test with p value 0.01
  if (type == "pptest" || type == "all")
  { 
    v=c(rep(0, p))
    for (h in p:1)
    {
      v[h] <- PP.test(Z[ , h], lshort = FALSE)$p.value
      if(v[h] > alpha) break
    }
    r_hat2 = p-h
    
    # older way
    # v=c(rep(0, p))
    # for (h in 1:p)
    # {
    #   v[h]=PP.test(Z[ , h],lshort = FALSE)$p.value
    # }
    # r_hat2 = p-sum(v > 0.01)
    
    
    
    if (type == "pptest"){
      METHOD <- c("Identifying cointegration rank of given time series",
                  "using Phillips-Perron test")
      names(r_hat2) <- "The estimated number of cointegration rank"
      names(lag.k) <- "Time lag"
      return(structure(list(Z = Z, coint_rank = r_hat2, lag.k = lag.k,
                            method = METHOD),
                       class = "coint"))
    }
  }
  # (d) Using unit-root test based on sample covariance with p value 0.01
  if (type == "Chang" || type == "all")
  {
    # older way
    # v=c(rep(0, p))
    # for (h in 1:p)
    # {
    #   v[h]=ur.test(Z[,h], lagk.vec=2, alpha=0.01)$result[1,1]
    # }
    # r_hat3 = p-sum(v==1)
    # print(v)
    # print(r_hat3)
    
    v=c(rep(0, p))
    for (l in p:1)
    {
      v[l] <- UR_test(Z[ , l], lagk.vec = 2, alpha = alpha)$reject[1, 1]
      if(v[l] == 1) {break}
    }
    r_hat3 = p - l
    if (type == "Chang"){
      METHOD <- c("Identifying cointegration rank of given time series",
                  "using Chang's method")
      names(r_hat3) <- "The estimated number of cointegration rank"
      names(lag.k) <- "Time lag"
      return(structure(list(Z = Z, coint_rank = r_hat3, lag.k = lag.k,
                            method = METHOD),
                       class = "coint"))
    }
  }
  if(type == "all"){
    result <- t(as.matrix(c(r_hat1,r_hat2,r_hat3)))
    rownames(result) <- "r_hat"
    colnames(result) <- c("acf","pptest","chang")
    METHOD <- c("Identifying cointegration rank of given time series",
                "using both three methods")
    names(lag.k) <-"Time lag"
    return(structure(list(Z = Z, coint_rank = result, lag.k = lag.k,
                          method = METHOD),
                     class = "coint"))
  }
}
