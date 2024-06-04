#' @name MartG_test
#' @title Testing for martingale difference hypothesis in high dimension
#' @description \code{MartG_test()} implements a new test proposed in
#'  Chang, Jiang and Shao (2021) for the following hypothesis testing problem: 
#' \deqn{H_0:\{{\bf x}_t\}_{t=1}^n\mathrm{\ is\ a\ MDS\ \ versus\ \ }H_1:
#' \{{\bf x}_t\}_{t=1}^n\mathrm{\ is\ not\ a\ MDS,} } where 
#' MDS is the abbreviation of "martingale difference sequence".
#' 
#' @param X \eqn{{\bf X} = \{{\bf x}_1, \dots , {\bf x}_n \}'}, an \eqn{n\times
#'   p} sample matrix, where \eqn{n} is the sample size and \eqn{p} is the 
#'   dimension of \eqn{{\bf x}_t}.
#' @param lag.k Time lag \eqn{K}, a positive integer, used to calculate the test
#'   statistic. Default is \code{lag.k} \eqn{=2}.
#' @param B Bootstrap times for generating multivariate normal distributed 
#' random vectors in calculating the critical value. 
#' Default is \code{B} \eqn{=1000}.
#' @param type String, a map is chosen by the \proglang{R} users, such as the
#'   default option is \code{'Linear'} means linear identity 
#'   map (\eqn{\boldsymbol \phi({\bf x})={\bf x}}). Also including another 
#'   option \code{'Quad'} (Both linear and quadratic terms 
#'   \eqn{\boldsymbol \phi({\bf x})=\{{\bf x}',({\bf x}^2)'\}'}). Also the users
#'   can choose set the map themselves, use for example \code{expression(X, X^2)},
#'    \code{quote(X, X^2)}, \code{parse(X, X^2)}, \code{substitute(X, X^2)} or 
#'    just map without function (such as cbind(X, X^2)) to set their own map. 
#'   See Section 2.1 in Chang, Jiang and Shao (2021) for more information.
#' @param alpha The prescribed significance level. Default is 0.05.
#' @param kernel.type String, an option for choosing the symmetric kernel 
#'                    used in the estimation of long-run covariance matrix, 
#'                    for example, \code{'QS'} (Quadratic spectral kernel), 
#'                    \code{'Par'} (Parzen kernel) and \code{'Bart'} 
#'                    (Bartlett kernel), see Andrews (1991) for more 
#'                    information. Default option is \code{kernel.type = 'QS'}.

#' @return An object of class "hdtstest" is a list containing the following
#'   components:
#'    
#'.  \item{statistic}{The value of the test statistic.}
#'   \item{p.value}{Numerical value which represents the p-value of the test.}
#'   \item{lag.k}{The time lag used in function.}
#'   \item{method}{A character string indicating what method was performed.}
#'   \item{type}{A character string which map used on data matrix \code{X}.}
#'   \item{kernel.type}{A character string indicating what kenel method was performed.}
#'   
#' @references Chang, J., Jiang, Q. & Shao, X. (2022). \emph{Testing the
#'   martingale difference hypothesis in high dimension}.  Journal of 
#'   Econometrics, in press
#' @examples
#' n <- 200
#' p <- 10
#' X <- matrix(rnorm(n*p),n,p)
#' res <- MartG_test(X, type="Linear")
#' res <- MartG_test(X, type=cbind(X, X^2)) #the same as Linear type
#' res <- MartG_test(X, type=quote(cbind(X, X^2))) # expr using quote
#' res <- MartG_test(X, type=substitute(cbind(X, X^2))) # expr using substitute
#' res <- MartG_test(X, type=expression(cbind(X, X^2))) # expr using expression
#' res <- MartG_test(X, type=parse(text="cbind(X, X^2)")) # expr using parse
#' map_fun <- function(X) {X <- cbind(X,X^2); X}
#' res <- MartG_test(X, type=map_fun)
#' Pvalue <- res$p.value
#' rej <- res$reject
#' @useDynLib HDTSA
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @importFrom methods formalArgs
#' @import Rcpp
#' @export


MartG_test <- function (X, lag.k=2, B=1000, type=c('Linear','Quad'), 
                      alpha=0.05, kernel.type=c('QS','Par','Bart')) {
  data_name <- all.vars(substitute(X))
  if (is.character(type)) {
    type <- match.arg(type)
    if (type == 'Linear') {
      Xj <- X
      d <- ncol(Xj)
      } 
    if (type == 'Quad'){
      Xj <- cbind(X,X^2)
      d <- ncol(Xj)
      }
    storage.mode(d) <- "integer"
  }
  else if (is.name(type)){
    varnames_expr <- all.vars(type)
    if (data_name != varnames_expr) 
      stop("expr is not validate, The variable name is different from the data name.")
    # print(type)
    Xj <- eval(type) 
    d <- ncol(Xj)
  }
  else if (is.call(type) || is.expression(type)){
    varnames_expr <- all.vars(type)
    if (length(varnames_expr) != 1)
      stop("expr is not validate, only support one data name and expr must include data name.")
    if (data_name != varnames_expr)
      stop("expr is not validate, The variable name is different from the data name.")
    #print(type)
    Xj <- eval(type)
    d <- ncol(Xj)
  }
  else if(is.function(type)){
    f <- type
    f_para <- formalArgs(f)
    varnames_expr <- f_para[1]
    if (data_name != varnames_expr)
      stop("function is not validate, The first argument of function is different
           from the data name.")
    if(!is.null(args) && length(f_para)>1)
      Xj <- do.call(f, args=c(list(X),args))
    else Xj <- f(X)
    d <- ncol(Xj)
  }
  else if(is.matrix(type)){
    Xj <- type
    names(type) <- "User define"
    d <- ncol(Xj)
  }
  else {
    expr <- substitute(type)
    varnames_expr <- all.vars(expr)
    if (length(varnames_expr) != 1)
      stop("expr is not validate, only support one data name.")
    if (data_name != varnames_expr)
      stop("expr is not validate, The variable name is different from the data name.")
   # print(expr)
    Xj <- eval(expr)
    d <- ncol(Xj)
  }
  
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
  
  
  # ---------- step 2: compute test statistic ----------
  Tn <- MartG_TestStatC(n, lag.k, X, Xj)
  
  ft <- MartG_ftC(n, lag.k, p, d, X, Xj)
  bn <- bandwith(ft,lag.k,p,d,ken_type)
	
	Gnstar <- MartG_bootc(n, lag.k, p, d, B, bn, ken_type, ft)  # critical value
	p.value <- mean(Gnstar>Tn)
	names(Tn) <- "Statistic"
	names(lag.k) <-"Time lag"
	names(kernel.type) <- "Symmetric kernel"
	METHOD <- "Testing for martingale difference hypothesis in high dimension"
	structure(list(statistic = Tn, p.value = p.value, lag.k=lag.k,
	               method = METHOD,
	               type = type, kernel = kernel.type),
	          class = "hdtstest")
}




