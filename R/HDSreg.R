#' @name HDSReg
#' @title High dimensional stochastic regression with latent factors
#' @description
#' \code{HDSReg()} consider a multivariate time series model which represents a high dimensional vector process as a sum of three terms: a linear 
#' regression of some observed regressors, a linear combination of some latent and serially correlated factors, and a vector white noise. \code{HDSReg} estimate the unknown 
#' regression coefficient matrix \eqn{D}, the number of factors \eqn{r} and the factor loading matrix \eqn{A}
#' @param Y A data matrix used for regression and factor inference with \eqn{n} rows and \eqn{p} columns, where \eqn{n} is the sample size and \eqn{p} is the dimension of the time series.
#' @param Z A data matrix with \eqn{n} rows and \eqn{m} columns, where \eqn{n} is the sample size and \eqn{m} is the dimension of the time series. it represents some observed regressors.
#' @param D A \eqn{p\times m} data matrix means the regression coefficient matrix \eqn{D}. Default is \code{NULL}, if \code{D} is not given, then the program will estimate it.
#' @param lag.k A time lag which is a positive integer specified to calculate \eqn{ \hat{\mathbf{M}}}.
#' @param twostep Logical. If \code{FALSE} (the default), then standard procedures for factor number inference (see \code{\link{factors}}) will be implemented. If \code{TRUE}, then a two step estimation procedure will be implemented 
#' for estimation of weak factors.
#' @seealso \code{\link{factors}}.
#' @return An object of class "HDSReg" is a list containing the following components:
#'
#' \item{factor_num}{Number of factors \eqn{r}.}
#' \item{reg.coff.mat}{The estimated \eqn{p \times m} regression coefficient matrix if \code{D} is not given.}
#' \item{loading.mat}{The estimated \eqn{p \times m} factor loading matrix.}
#' @references Chang, J., Guo, B. & Yao, Q. (2015).  \emph{High dimensional stochastic regression with latent factors, endogeneity and nonlinearity}, Journal of Econometrics, Vol. 189, pp. 297–312.
#' @examples 
#' n <- 400
#' p <- 200
#' m <- 2
#' r <- 3
#' X <- mat.or.vec(n,r)
#' x1 <- arima.sim(model=list(ar=c(0.6)),n=n)
#' x2 <- arima.sim(model=list(ar=c(-0.5)),n=n)
#' x3 <- arima.sim(model=list(ar=c(0.3)),n=n)
#' X <- cbind(x1,x2,x3)
#' X <- t(X)
#' 
#' Z <- mat.or.vec(m,n)
#' S1 <- matrix(c(5/8,1/8,1/8,5/8),2,2)
#' Z[,1] <- c(rnorm(m))
#' for(i in c(2:n)){
#'   Z[,i] <- S1%*%Z[, i-1] + c(rnorm(m))
#' }
#' D <- matrix(runif(p*m, -2, 2), ncol=m)
#' A <- matrix(runif(p*r, -2, 2), ncol=r)
#' eps <- mat.or.vec(n, p)
#' eps <- matrix(rnorm(n*p), p, n)
#' Y <- D %*% Z + A %*% X + eps
#' Y <- t(Y)
#' Z <- t(Z)
#' res1 <- HDSReg(Y,Z,D,lag.k=2)
#' res2 <- HDSReg(Y,Z,lag.k=2)
#' @useDynLib HDTSA
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @import Rcpp 
#' @import RcppEigen
#' @export
HDSReg <- function (Y,Z,D=NULL,lag.k=1,twostep=FALSE) {
  Y_n <- nrow(Y)
  Y_p <- ncol(Y)
  Z_n <- nrow(Z)
  Z_p <- ncol(Z)
  storage.mode(Y_n) <- "integer"
  storage.mode(Y_p) <- "integer"
  storage.mode(Z_n) <- "integer"
  storage.mode(Z_p) <- "integer"
  
  #first step: estimate D (Least square) #without endogeneity and nolinearity
  if(missing(D)) {
    D <- solve(MatMult(t(Z), Z))
    D <- MatMult(D, MatMult(t(Z), Y))
    D <- t(D)
  }
  
  #estimate yt=Dzt+nt
  
  eta <- t(Y)-MatMult(D, t(Z))
  eta <- t(eta)
  factor_list <- factors(eta, lag.k, twostep)
  r <- factor_list$factor_num
  loading.mat <- factor_list$loading.mat
  outlist <- list(reg.coff.mat=D, factor_num=r, loading.mat=loading.mat)
  class(outlist) <- c("HDSReg")
  return(outlist)
  
}
