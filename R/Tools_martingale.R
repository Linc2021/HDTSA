#' @name MartG_test
#' @title Testing for martingale difference hypothesis in high dimension
#' @description \code{MartG_test()} implements a new test proposed in
#'  Chang, Jiang and Shao (2023) for the following hypothesis testing problem: 
#' \deqn{H_0:\{{\bf y}_t\}_{t=1}^n\mathrm{\ is\ a\ MDS\ \ versus\ \ }H_1:
#' \{{\bf y}_t\}_{t=1}^n\mathrm{\ is\ not\ a\ MDS}\,, } where 
#' MDS is the abbreviation of "martingale difference sequence".
#' 
#' @details
#' Write \eqn{{\bf x}= (x_1,\ldots,x_p)'}.
#' When \code{type = "Linear"}, the linear identity map is defined
#' as \eqn{\boldsymbol \phi({\bf x})={\bf x}}.
#' 
#' When \code{type = "Quad"},
#' \eqn{\boldsymbol \phi({\bf x})=\{{\bf x}',({\bf x}^2)'\}'}
#' includes both linear and quadratic terms, where
#' \eqn{{\bf x}^2 = (x_1^2,\ldots,x_p^2)'}.
#' 
#' We can also choose \eqn{\boldsymbol \phi({\bf x}) = \cos({\bf x})} to capture
#' certain type of nonlinear dependence, where
#' \eqn{\cos({\bf x}) = (\cos x_1,\ldots,\cos x_p)'}.
#' 
#' See 'Examples'.
#' 
#' 
#' 
#' 
#' @param Y An \eqn{n \times p} data matrix \eqn{{\bf Y} = ({\bf y}_1, \dots , {\bf y}_n )'},
#'  where \eqn{n} is the number of the observations of the \eqn{p \times 1}
#'  time series \eqn{\{{\bf y}_t\}_{t=1}^n}.
#' @param lag.k The time lag \eqn{K} used to calculate the test
#'   statistic [See (3) in Chang, Jiang and Shao (2023)]. The default is 2.
#' @param B The number of bootstrap replications for generating multivariate
#' normally distributed random vectors when calculating the critical value. 
#' The default is 1000. 
#' @param type The map used for constructing the test statistic.
#' Available options include: \code{"Linear"} (the default) for the linear identity 
#' map and \code{"Quad"} for the map including both linear and quadratic terms.
#' \code{type} can also be set by the users. See 'Details' and Section 2.1
#' of Chang, Jiang and Shao (2023) for more information.
#' @param alpha The significance level of the test. The default is 0.05.
#' @param kernel.type The option for choosing the symmetric kernel 
#'                    used in the estimation of long-run covariance matrix. 
#'                    Available options include: \code{"QS"} (the default)  for the
#'                    Quadratic spectral kernel, \code{"Par"} for the Parzen kernel,
#'                    and \code{"Bart"}  for the Bartlett kernel.
#'                    See Chang, Jiang and Shao (2023) for more information.

#' @return An object of class \code{"hdtstest"}, which contains the following
#'   components:
#'   \item{statistic}{The test statistic of the test.}
#'   \item{p.value}{The p-value of the test.}
#'   \item{lag.k}{The time lag used in function.}
#'   \item{type}{The map used in function.}
#'   \item{kernel.type}{The kernel used in function.}
#'   
#' @references Chang, J., Jiang, Q., & Shao, X. (2023). Testing the
#'   martingale difference hypothesis in high dimension. \emph{Journal of 
#'   Econometrics}, \strong{235}, 972--1000. \doi{doi:10.1016/j.jeconom.2022.09.001}.
#' @examples
#' # Example 1
#' n <- 200
#' p <- 10
#' X <- matrix(rnorm(n*p),n,p)
#' 
#' res <- MartG_test(X, type="Linear")
#' res <- MartG_test(X, type=cbind(X, X^2)) #the same as type = "Quad"
#' 
#' ## map can also be defined as an expression in R.
#' res <- MartG_test(X, type=quote(cbind(X, X^2))) # expr using quote()
#' res <- MartG_test(X, type=substitute(cbind(X, X^2))) # expr using substitute()
#' res <- MartG_test(X, type=expression(cbind(X, X^2))) # expr using expression()
#' res <- MartG_test(X, type=parse(text="cbind(X, X^2)")) # expr using parse()
#' 
#' ## map can also be defined as a function in R.
#' map_fun <- function(X) {X <- cbind(X, X^2); X}
#' 
#' res <- MartG_test(X, type=map_fun)
#' Pvalue <- res$p.value
#' rej <- res$reject
#' @useDynLib HDTSA
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @importFrom methods formalArgs
#' @export


MartG_test <- function (Y, lag.k = 2, B = 1000, type = c("Linear", "Quad"), 
                      alpha = 0.05, kernel.type = c("QS", "Par", "Bart")) {
  data_name <- all.vars(substitute(Y))
  if (is.character(type)) {
    type <- match.arg(type)
    if (type == 'Linear') {
      Yj <- Y
      d <- ncol(Yj)
      } 
    if (type == 'Quad'){
      Yj <- cbind(Y, Y^2)
      d <- ncol(Yj)
      }
  }
  else if (is.name(type)){
    varnames_expr <- all.vars(type)
    if (data_name != varnames_expr) 
      stop("expr is not validate, The variable name is different from the data name.")
    # print(type)
    Yj <- eval(type) 
    d <- ncol(Yj)
    type <- deparse(type)
  }
  else if (is.call(type) || is.expression(type)){
    varnames_expr <- all.vars(type)
    if (length(varnames_expr) != 1)
      stop("expr is not validate, only support one data name and expr must include data name.")
    if (data_name != varnames_expr)
      stop("expr is not validate, The variable name is different from the data name.")
    #print(type)
    Yj <- eval(type)
    d <- ncol(Yj)
    if(is.call(type)) type <- deparse(type)
    if(is.expression(type)) type<- deparse(type[[1]])
  }
  else if(is.function(type)){
    f <- type
    f_para <- formalArgs(f)
    varnames_expr <- f_para[1]
    if (data_name != varnames_expr)
      stop("function is not validate, The first argument of function is different
           from the data name.")
    if(!is.null(args) && length(f_para)>1)
      Yj <- do.call(f, args=c(list(Y),args))
    else Yj <- f(Y)
    d <- ncol(Yj)
    type <- "User-defined"
  }
  else if(is.matrix(type)){
    Yj <- type
    # names(type) <- "User define"
    d <- ncol(Yj)
    type <- "User-defined"
  }
  else {
    expr <- substitute(type)
    varnames_expr <- all.vars(expr)
    if (length(varnames_expr) != 1)
      stop("expr is not validate, only support one data name.")
    if (data_name != varnames_expr)
      stop("expr is not validate, The variable name is different from the data name.")
   # print(expr)
    Yj <- eval(expr)
    d <- ncol(Yj)
    type <- deparse(type)
  }
  
  kernel.type <- match.arg(kernel.type)
  ken_type <- switch(kernel.type,
                        "QS" = 1,
                        "Par" = 2,
                        "Bart" = 3)
  
  n <- nrow(Y)
  p <- ncol(Y)
  
  # ---------- step 1: Transformation function ----------
  
  
  # ---------- step 2: compute test statistic ----------
  Tn <- MartG_TestStatC(n, lag.k, Y, Yj)
  
  ft <- MartG_ftC(n, lag.k, p, d, Y, Yj)
  bn <- bandwith(ft,lag.k,p,d,ken_type)
  
  boot_nomal <- matrix(rnorm(B*(n-lag.k)), B, n-lag.k)
  Gnstar <- MartG_bootc(n, lag.k, p, d, B, bn, ken_type, ft, boot_nomal)  # critical value
	p.value <- mean(Gnstar>Tn)
	names(Tn) <- "Statistic"
	names(lag.k) <-"Time lag"
	names(kernel.type) <- "Symmetric kernel"
	structure(list(statistic = Tn, p.value = p.value, lag.k=lag.k,
	               type = type, kernel = kernel.type),
	          class = "hdtstest")
}




