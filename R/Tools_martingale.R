#' @name MartG_test
#' @title Testing for martingale difference hypothesis in high dimension
#' @description \code{MartG_test()} implements a new test for the following martingale
#' difference hypothesis for high-dimensional time series: \deqn{H_0:\{{\bf x}_t
#' \}_{t=1}^n\mathrm{\ is\ a\ MDS\ \ versus\ \ }H_1:\{{\bf x}_t
#' \}_{t=1}^n\mathrm{\ is\ not\ a\ MDS,} } where MDS is the abbreviation of "martingale difference
#' sequence". The test (Chang, Jiang and Shao, 2021) is built on the sup-norm of a matrix-valued sample
#' nonlinear dependencemetric at a finite and possibly growing number of lags.
#' To the best of knowledge, this is the first valid test for the martingale
#' difference hypothesis that allows for large dimension and captures nonlinear
#' serial dependence. Besides, the test is robust to conditional moments of
#' unknown forms and panel dependence of unknown magnitude.
#' 
#' @param X \eqn{{\bf X} = \{{\bf x}_1, \dots , {\bf x}_n \}'}, an \eqn{n\times
#'   p} sample matrix, where \eqn{n} is the sample size and \eqn{p} is the dimension of \eqn{{\bf x}_t}.
#' @param lag.k Time lag \eqn{K}, a positive integer, used to calculate test
#'   statistic. Default is \code{lag.k} \eqn{=2}.
#' @param B Bootstrap times for generating vector from a multivariate normal
#'   distribution. Default is \code{B} \eqn{=2000}.
#' @param type String, a map is chosen by the \proglang{R} users, such as the
#'   default option is \code{'Linear'} means linear identity map (\eqn{\boldsymbol \phi({\bf x})={\bf x}}). Also including
#'   another option \code{'Quad'} (Both linear and quadratic terms \eqn{\boldsymbol \phi({\bf x})=\{{\bf x}',({\bf x}^2)'\}'}).
#' @param alpha The significance level used for testing. Default is 0.05.
#' @param kernel.type String, an option for choosing an optimal symmetric kernel
#'   type, for example, \code{'QS'} (Quadratic Spectral), \code{'Par'} (Parzen), \code{'Bart'} (Bartlett), see Andrews (1991) for more information. Default option is
#'   \code{kernel.type = 'QS'}.

#' @return An object of class "MartG_test" is a list containing the following components:
#'
#' \item{reject}{Logical value. If \code{TRUE}, it means rejecting the null hypothesis, otherwise it means not rejecting the null hypothesis }
#' \item{p.value}{Numerical value which represents the p-value of the test based on the observed data \eqn{\{{\bf x}_t\}_{t=1}^n}}
#' @references Chang, J., Jiang, Q. & Shao, X. (2021). \emph{Testing the martingale difference hypothesis in high dimension}.
#' @examples 
#' n <- 200
#' p <- 150
#' X <- matrix(rnorm(n*p),n,p)
#' res <- MartG_test(X)
#' Pvalue <- res$p.value
#' rej <- res$reject
#' @useDynLib HDTSA
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @import Rcpp 
#' @import RcppEigen
#' @export


MartG_test <- function (X, lag.k=2, B=1000, type=c('Linear','Quad','In'), 
                      alpha=0.05, kernel.type=c('QS','Par','Bart')) {
  
  type <- match.arg(type)
  kernel.type <- match.arg(kernel.type)
  ken_type <- switch(kernel.type,
                        "QS" = 1,
                        "Par" = 2,
                        "Bart" = 3)
  
  n <- nrow(X)
  p <- ncol(X)
  storage.mode(p) <- "integer"
  storage.mode(n) <- "integer"
  
  # ---------- step 1: Transformation function ----------
  if (type == 'Linear') {
    Xj <- X
    d <- ncol(Xj)
  } 
  if (type == 'Quad'){
    Xj <- cbind(X,X^2)
    d <- ncol(Xj)
  }
  if (type == 'In'){
    InX <- matrix(NA,nrow = n, ncol = p*(p-1)/2)
    flag <- 1
    for (i  in c(1:(p-1))) {
      for (j in c((i+1):p)) {
        InX[,flag] <- X[,i] * X[,j]
        flag <- flag + 1
      }
    }
    Xj <- cbind(X, X^2, InX)
    d <- ncol(Xj)
  }
  storage.mode(d) <- "integer"
  
  # ---------- step 2: compute test statistic ----------
  Tn <- MartG_TestStatC(n, lag.k, X, Xj)
  
  ft <- MartG_ftC(n, lag.k, p, d, X, Xj)
  bn <- MartG_bandwith(ft,lag.k,p,d,ken_type)
	
	Gnstar <- MartG_bootc(n, lag.k, p, d, B, bn, ken_type, ft)  # critical value
	p.value <- mean(Gnstar>Tn)
	Results <- list(reject = (p.value<0.05), p.value = p.value)
	class(Results) <- c("MartG_test")
	return(Results)
}




