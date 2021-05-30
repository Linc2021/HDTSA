#' @name factors
#' @title Factor modeling: Inference for the number of factors
#' @description
#' \code{factors()} deals with factor modeling for high-dimensional time series based on a dimension-reduction 
#' viewpoint proposed in Lam and Yao (2012):\deqn{{\bf y}_t = {\bf Ax}_t + {\boldsymbol{\epsilon}}_t, } where \eqn{{\bf x}_t} is an \eqn{r \times 1} latent process with (unknow) \eqn{r \leq p},
#' \eqn{{\bf A}} is a \eqn{p \times r} unknown constant matrix, and \eqn{ {\boldsymbol{\epsilon}}_t \sim \mathrm{WN}({\boldsymbol{\mu}}_{\epsilon}, {\bf \Sigma}_{\epsilon})} is a vector white-noise process. 
#' Under stationary settings, the inference is simple in the sense that both the number of factors \eqn{r} and the factor loadings \eqn{{\bf A}} can be estimated in terms of an eigenanalysis for a nonnegative definite matrix, and is 
#' therefore applicable when the dimension of time series is on the order of a few thousands.
#' 
#' @param Y \eqn{{\bf Y} = \{{\bf y}_1, \dots , {\bf y}_n \}'}, a data matrix used for factor inference with \eqn{n} rows and \eqn{p} columns, where \eqn{n} is the sample size and \eqn{p} is the dimension of the time series.
#' @param lag.k Time lag \eqn{k_0} used to calculate the nonnegative definte matrix \eqn{ \widehat{\mathbf{M}}}: \deqn{\widehat{\mathbf{M}}\ =\ \sum_{k=1}^{k_0}\widehat{\mathbf{\Sigma}}_y(k)\widehat{\mathbf{\Sigma}}_y(k)', }
#'              where \eqn{\widehat{\bf \Sigma}_y(k)} is the sample autocovariance of \eqn{ {\bf y}_t} at lag \eqn{k}.
#' @param twostep logical. If \code{FALSE} (the default), then standard procedures for factor number inference will be implemented. If \code{TRUE}, then a two step estimation procedure will be implemented 
#' for the estimation of weak factors.

#' @return An object of class "factors" is a list containing the following components:
#' \item{factor_num}{The estimated Number of factors \eqn{\hat{r}}}
#' \item{loading.mat}{The estimated \eqn{p \times r} factor loading matrix \eqn{\widehat{\bf A}}}
#' @references Lam, C. & Yao, Q. (2012). \emph{Factor modelling for high-dimensional time series: Inference for the number of factors}. The Annals of Statistics, 40, 694-726.
#' @examples
#' ## Generate x_t 
#' p <- 400
#' n <- 400
#' r <- 3
#' X <- mat.or.vec(n, r)
#' A <- matrix(runif(p*r, -1, 1), ncol=r)
#' x1 <- arima.sim(model=list(ar=c(0.6)), n=n)
#' x2 <- arima.sim(model=list(ar=c(-0.5)), n=n)
#' x3 <- arima.sim(model=list(ar=c(0.3)), n=n)
#' eps <- matrix(rnorm(n*p), p, n)
#' X <- t(cbind(x1, x2, x3))
#' Y <- A %*% X + eps
#' Y <- t(Y)
#' fac <- factors(Y,lag.k=2)
#' r_hat <- fac$factor_num
#' loading_Mat <- fac$loading.mat
#' @useDynLib HDTSA
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @import Rcpp 
#' @import RcppEigen
#' @export



factors <- function (Y,lag.k=5,twostep=FALSE) {
  
  n <- nrow(Y)
  p <- ncol(Y)
  r <- 0
  step <- ifelse(twostep==TRUE,2,1)
  Y <- t(Y)
  
  storage.mode(r) <- "integer"
  storage.mode(p) <- "integer"
  storage.mode(n) <- "integer"
  storage.mode(step) <- "integer"
  storage.mode(Y) <- "double"
  #step1 caculate the M
  
  for(time in c(1:step)){
    #two method while estimate M
    
    M <- diag(rep(0,p))
    mean_y <- as.matrix(rowMeans(Y))
    
    storage.mode(M) <- "double"
    storage.mode(mean_y) <- "double"
    
    for(k in 1:lag.k){
      Sigma_y <- sigmak(Y, mean_y, k, n)
      M <- M + MatMult(Sigma_y, t(Sigma_y))
    }
    
    #step2 Eigendecoposition
    t <- eigen(M, symmetric=T)
    ev <- t$value
    G <- as.matrix(t$vectors)
    
    p1 <- ceiling(p/4)
    ratio <- ev[2:(p1+1)] / ev[1:p1]
    min_ratio <- min(ratio)
    index <- 1:p1
    r <- r + index[ratio==min_ratio]
    if(time == 1){
      final_vector <- G[,c(1:r)]
      if(twostep){
        Y <- Y - final_vector %*% t(final_vector) %*% Y
      }
    }
    else{
      final_vector <- cbind(final_vector, G[, c(1:r)])
    }
  }
  outlist <- list(factor_num = r, loading.mat = final_vector)
  class(outlist) <- c("factors")
  return(outlist)
    
    #extension: two step method
    

}
