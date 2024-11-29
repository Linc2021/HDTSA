#' @name Factors
#' @title Factor analysis for vector time series
#' @description \code{Factors()} deals with factor modeling for high-dimensional
#' time series proposed in Lam and Yao (2012):\deqn{{\bf y}_t = {\bf Ax}_t +
#' {\boldsymbol{\epsilon}}_t, } where \eqn{{\bf x}_t} is an \eqn{r \times 1}
#' latent process with (unknown) \eqn{r \leq p}, \eqn{{\bf A}} is a \eqn{p
#' \times r} unknown constant matrix, and \eqn{ {\boldsymbol{\epsilon}}_t} is a
#' vector white noise process. The number of factors \eqn{r} and the factor
#' loadings \eqn{{\bf A}} can be estimated in terms of an eigenanalysis for a
#' nonnegative definite matrix, and is therefore applicable when the dimension
#' of \eqn{{\bf y}_t} is on the order of a few thousands. This function aims to
#' estimate the number of factors \eqn{r} and the factor loading matrix
#' \eqn{{\bf A}}.
#' 
#' @details
#' The threshold operator \eqn{T_\delta(\cdot)} is defined as
#' \eqn{T_\delta({\bf W}) = \{w_{i,j}1(|w_{i,j}|\geq \delta)\}} for any matrix
#' \eqn{{\bf W}=(w_{i,j})}, with the threshold level \eqn{\delta \geq 0} and \eqn{1(\cdot)}
#' representing the indicator function. We recommend to choose
#'   \eqn{\delta=0} when \eqn{p} is fixed and \eqn{\delta>0} when \eqn{p \gg n}.
#'
#' @param Y An \eqn{n \times p} data matrix \eqn{{\bf Y} = ({\bf y}_1, \dots , {\bf y}_n )'},
#' where \eqn{n} is the number of the observations of the \eqn{p \times 1} time
#' series \eqn{\{{\bf y}_t\}_{t=1}^n}.
#' @param lag.k The time lag \eqn{K} used to calculate the nonnegative definite
#'   matrix \eqn{ \hat{\mathbf{M}}}: \deqn{\hat{\mathbf{M}}\ =\
#'   \sum_{k=1}^{K} T_\delta\{\hat{\mathbf{\Sigma}}_y(k)\} T_\delta\{\hat{\mathbf{\Sigma}}_y(k)\}'\,,
#'   } where \eqn{\hat{\bf \Sigma}_y(k)} is the sample autocovariance of
#'   \eqn{ {\bf y}_t} at lag \eqn{k} and \eqn{T_\delta(\cdot)}
#'   is a threshold operator with the threshold level \eqn{\delta \geq 0}. See 'Details'.
#'   The default is 5.
#' @param twostep Logical. If \code{twostep = FALSE} (the default), the standard
#'   procedure [See Section 2.2 in Lam and Yao (2012)] for estimating \eqn{r}
#'   and \eqn{{\bf A}} will be implemented. If \code{twostep = TRUE}, the two-step
#'   estimation procedure [See Section 4 in Lam and Yao (2012)]
#'   for estimating \eqn{r} and \eqn{{\bf A}} will be implemented.
#' @param thresh  Logical. If \code{thresh = FALSE} (the default), no thresholding will
#'   be applied to estimate \eqn{\hat{\mathbf{M}}}. If \code{thresh = TRUE},
#'   \eqn{\delta} will be set through \code{delta}.
#' @param delta  The value of the threshold level \eqn{\delta}. The default is
#'  \eqn{ \delta = 2 \sqrt{n^{-1}\log p}}.

#' @return An object of class \code{"factors"}, which contains the following
#'   components: 
#'   \item{factor_num}{The estimated number of factors
#'   \eqn{\hat{r}}.} 
#'   \item{loading.mat}{The estimated \eqn{p \times \hat{r}} factor
#'   loading matrix \eqn{\hat{\bf A}}.}
#'   \item{X}{The \eqn{n\times \hat{r}} matrix
#'   \eqn{\hat{\bf X}=(\hat{\bf x}_1,\dots,\hat{\bf x}_n)'} with
#'   \eqn{\hat{\bf x}_t = \hat{\bf A}'\hat{\bf y}_t}.}
#'   \item{lag.k}{The time lag used in function.}
#'   
#'  

#' @references Lam, C., & Yao, Q. (2012). Factor modelling for
#'   high-dimensional time series: Inference for the number of factors. \emph{The
#'   Annals of Statistics}, \strong{40}, 694--726. \doi{doi:10.1214/12-AOS970}.
#' @examples
#' # Example 1 (Example in Section 3.3 of lam and Yao 2012)
#' ## Generate y_t
#' p <- 200
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
#' 
#' fac <- Factors(Y,lag.k=2)
#' r_hat <- fac$factor_num
#' loading_Mat <- fac$loading.mat
#' @useDynLib HDTSA
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @export



Factors <- function (Y, lag.k = 5, thresh = FALSE, delta = 2 * sqrt(log(ncol(Y)) / nrow(Y)),
                     twostep = FALSE) 
{
  n <- nrow(Y)
  p <- ncol(Y)
  r <- c()
  Y1 <- Y
  step <- ifelse(twostep==TRUE, 2, 1)
  Y <- t(Y)
  
  #step1 caculate the M
  for(time in c(1:step)){
    #two method while estimate M
    
    M <- diag(rep(0,p))
    mean_y <- as.matrix(rowMeans(Y))
    
    
    for(k in 1:lag.k){
      Sigma_y <- sigmak(Y, mean_y, k, n)
      if(!thresh)
        M <- M + MatMult(Sigma_y, t(Sigma_y))
      else {
        Sigma_y <- thresh_C(Sigma_y, delta)
        M <- M + MatMult(Sigma_y, t(Sigma_y))
      }
      
    }
    
    #step2 Eigendecoposition
    t <- eigen(M, symmetric=T)
    ev <- t$value
    G <- as.matrix(t$vectors)
    
    p1 <- ceiling(p*0.75)
    ratio <- ev[2:(p1+1)] / ev[1:p1]
    min_ratio <- min(ratio)
    index <- 1:p1
    # r <- r + index[ratio==min_ratio]
    r <- c(r, which.min(ratio))
    if(time == 1){
      final_vector <- G[,c(1:r[time]), drop = FALSE]
      if(twostep){
        Y <- Y - MatMult(MatMult(final_vector, t(final_vector)), Y)
      }
    }
    else{
      final_vector <- cbind(final_vector, G[, c(1:r[time])])
    }
  }
  X <- MatMult(Y1, final_vector)
  
  # METHOD <- "Inference for the number of factors"
  names(lag.k) <-"Time lag"
  structure(list(factor_num = r, loading.mat = final_vector, X = X,
                 lag.k=lag.k),
            class = "factors")
}
