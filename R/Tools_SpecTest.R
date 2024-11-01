#' @name SpecTest
#' @title Statistical inference for high-dimensional spectral density matrix
#' @description \code{SpecTest()} implements a new global test proposed in 
#' Chang et al. (2023) for the following hypothesis testing problem: 
#' \deqn{H_0:f_{i,j}(\omega)=0 \mathrm{\ for\ any\ }(i,j)\in\mathcal{I}\mathrm{\ and\ }
#' \omega \in \mathcal{J}\mathrm{\ \ versus\ \ H_1:H_0\ is\ not\ true.} }
#' Here, \eqn{f_{i,j}(\omega)} represents the cross-spectral density between 
#' components \eqn{i} and \eqn{j} of a multivariate time series at frequency 
#' \eqn{\omega}. The set \eqn{\mathcal{I}} refers to the pairs of indices of 
#' interest (i.e., those under investigation for potential cross-spectral 
#' relationships), while \eqn{\mathcal{J}} represents the frequency domain over 
#' which the test is performed. The null hypothesis \eqn{H_0} posits that the 
#' cross-spectral density is zero for all specified pairs \eqn{(i,j)\in \mathcal{I}} and 
#' frequencies \eqn{\omega \in \mathcal{J}}, implying no linear relationship between these 
#' components across the frequency domain.
#' 
#' 
#' @param X \eqn{{\bf X} = \{{\bf x_1}, \dots , {\bf x}_n \}}, a \eqn{n\times
#'   p} sample matrix, where \eqn{n} is the sample size and \eqn{p} is the 
#'   dimension of \eqn{{\bf x}_t}.
#' @param B The bootstrap times for generating multivariate normal distributed 
#' random vectors in calculating the critical value. Default is \code{B} \eqn{=2000}.
#' @param flag_c The bandwidth \eqn{c} of the flat-top kernel for estimating 
#' \eqn{f_{i,j}(\omega)}, where \eqn{c\in(0,1]}. Default is \code{flag_c} \eqn{=0.8}.
#' @param J.set The set \eqn{\mathcal{J}} for frequencies, a vector, used to calculate the test statistic.
#' @param cross.indices The set \eqn{\mathcal{I}} for \eqn{(i,j)} is a matrix with two columns,
#'  employed in calculating the test statistic.
#' @seealso \code{\link{SpecMulTest}}
#' 
#' @return An object of class "hdtstest" is a list containing the following
#'   components:
#'
#'   \item{Stat}{Numerical value which represents the value of test statistic.}
#'   \item{pval}{Numerical value which represents the p-value of the test.}
#'   \item{cri95}{Numerical value which represents the critical value of the test
#'   at the significance level 0.05.}
#'   \item{method}{A character string indicating what method was performed.}
#' @references Chang, J., Jiang, Q., McElroy, T. S., & Shao, X. (2023). 
#' \emph{Statistical inference for high-dimensional spectral density matrix}. arXiv preprint arXiv:2212.13686.
#' 
#' @examples
#' n <- 200
#' p <- 10
#' flag_c <- 0.8
#' B <- 1000
#' burn <- 1000
#' z.sim <- matrix(rnorm((n+burn)*p),p,n+burn)
#' phi.mat <- 0.4*diag(p)
#' x.sim <- phi.mat %*% z.sim[,(burn+1):(burn+n)]
#' x <- x.sim - rowMeans(x.sim)
#' cross.indices <- matrix(c(1,2), ncol=2)
#' J.set <- 2*pi*seq(0,3)/4 - pi
#' res <- SpecTest(t(x), J.set, cross.indices, B, flag_c)
#' Stat <- res$statistic
#' Pvalue <- res$p.value
#' CriVal <- res$cri95
#' @export
SpecTest <- function(X, J.set, cross.indices, B = 1000, flag_c = 0.8)
{
  p <- ncol(X)
  n <- nrow(X)
  K <- length(J.set)
  r <- dim(cross.indices)[1]
  l.band <- ChooseLn(X, 2, 5) 
  X <- t(X)
  
  if(n > 2*l.band+1){
    ## ------------------------------------
    ##  Part II: computing Test distributions
    ## ------------------------------------
    GhatC <- CmpGammaC(X)  # List(n) : p*p 
    T2 <- TestStatC(GhatC,n,p,r,K,cross.indices,J.set,l.band,flag_c)  
    
    
    # ===================================================================
    #  critical value for studentized test statistic based on C.hat
    # ===================================================================
    
    n.tilde <- n - 2*l.band -1
    T2.stars = TestStarC(X, GhatC, n.tilde, n, p, r, K, flag_c, cross.indices, 
                         J.set, l.band, B, type=1)
    
    names(T2) <- "Statistic"
    cri95=sort(T2.stars)[floor(.95*B)]
    names(cri95) <- "the critical value of the test at the significance level 0.05"
    METHOD = "Statistical inference for high-dimensional spectral density matrix"
    structure(list(statistic=T2, p.value=sum(T2.stars>T2)/B, 
                 cri95=cri95,
                 method = METHOD),
                 class = "hdtstest")
  }else{
    warning("n<=2ln+1 which does not satisfy the requirement of this test.")
    structure(list(statistic = NA, p.value = NA, 
                   cri95 = NA,
                   method = METHOD),
              class = "hdtstest")
  }
}



#' @name SpecMulTest
#' @title Statistical inference for high-dimensional spectral density matrix
#' @description \code{SpecMulTest()} implements a new multiple test proposed in
#'  Chang et al. (2023) for the \eqn{Q} hypothesis testing problems: 
#' \deqn{H_{0,q}:f_{i,j}(\omega)=0\mathrm{\ for\ any\ }(i,j)\in\mathcal{I}^{(q)}\mathrm{\ and\ }
#' \omega\in\mathcal{J}^{(q)}\mathrm{\ \ versus\ \ }
#' H_{1,q}:H_{0,q}\mathrm{\ is\ not\ true.} }
#' for \eqn{q\in\{1,\dots,Q\}}. 
#' Here, \eqn{f_{i,j}(\omega)} represents the cross-spectral density between 
#' components \eqn{i} and \eqn{j} at frequency \eqn{\omega}. The sets 
#' \eqn{\mathcal{I}^{(q)}} and \eqn{\mathcal{J}^{(q)}} refer to the index pairs 
#' and frequency domain associated with the \eqn{q}-th test. The goal is to 
#' simultaneously test multiple hypotheses regarding the nullity of cross-spectral 
#' densities across different pairs and frequencies.
#' 
#' 
#' @param Q The number of the hypothesis tests. 
#' @param PVal A vector of length \code{Q} representing P-values for the \code{Q} hypothesis tests.
#' @param alpha The prescribed significance level. Default is 0.05.
#' @param seq_len Length used to take discrete points between 0 and 
#' \eqn{\sqrt{(2\times\log(Q)-2\times\log(\log(Q))}}. Default is 0.01.
#' @seealso \code{\link{SpecTest}}
#' 
#' @return An object of class "hdtstest" is a list containing the following
#'   components:
#'   \item{MultiTest}{Logical vector with length Q. If the element is \code{TRUE}, 
#'   it means rejecting the corresponding sub-null hypothesis, otherwise it means 
#'   not rejecting the corresponding sub-null hypothesis.}
#'   
#' @references Chang, J., Jiang, Q., McElroy, T. S., & Shao, X. (2023). 
#' \emph{Statistical inference for high-dimensional spectral density matrix}. arXiv preprint arXiv:2212.13686.
#' @examples
#' n <- 200
#' p <- 10
#' flag_c <- 0.8
#' B <- 1000
#' burn <- 1000
#' z.sim <- matrix(rnorm((n+burn)*p),p,n+burn)
#' phi.mat <- 0.4*diag(p)
#' x.sim <- phi.mat %*% z.sim[,(burn+1):(burn+n)]
#' x <- x.sim - rowMeans(x.sim)
#' Q <- 4
#' ISET <- list()
#' ISET[[1]] <- matrix(c(1,2),ncol=2)
#' ISET[[2]] <- matrix(c(1,3),ncol=2)
#' ISET[[3]] <- matrix(c(1,4),ncol=2)
#' ISET[[4]] <- matrix(c(1,5),ncol=2)
#' JSET <- as.list(2*pi*seq(0,3)/4 - pi)
#' PVal <- rep(NA,Q)
#' for (q in 1:Q) {
#'   cross.indices <- ISET[[q]]
#'   J.set <- JSET[[q]]
#'   temp.q <- SpecTest(t(x), J.set, cross.indices, B, flag_c)
#'   PVal[q] <- temp.q$p.value
#' }  # Q
#' res <- SpecMulTest(Q, PVal)
#' res
#' 
#' @useDynLib HDTSA
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @import Rcpp
#' @export
SpecMulTest <- function(Q, PVal, alpha=0.05, seq_len=0.01){
  
  Tseq = seq(0,sqrt(2*log(Q)-2*log(log(Q))),seq_len)
  
  Temp <- matrix(NA, Q, length(Tseq))
  for (kk in 1:length(Tseq)) {
    Temp[,kk] <- as.matrix(qnorm(1-PVal) >= Tseq[kk]) # 1 * Q
  }
  FDP <- Q*(1-pnorm(Tseq)) / pmax(apply(Temp, 2, sum), 1)
  if (length(which(FDP<=alpha)) > 0){
    hatT <- min(Tseq[which(FDP<=alpha)])
    t_FDR <- hatT
  }else{
    t_FDR <- sqrt(2*log(Q))
  }
  
  # compute FDR
  MultiTest <- (qnorm(1-PVal) >= t_FDR) # Q
  
  METHOD = "Statistical inference for high-dimensional spectral density matrix"
  
  structure(list(MultiTest=MultiTest,
                         method = METHOD),
                    class = "hdtstest")
}



# selection of ln
ChooseLn <- function(x, cn, Kn){
  n <- nrow(x)
  p <- ncol(x)
  gamma0 <- diag((diag(t(x) %*% x))^(-1/2),nrow=p)
  rhoK <- rep(NA, Kn)
  for (k in 2:(Kn+1)) {
    rhoK[k-1] <- mean( abs( gamma0 %*% (t(x[(k+1):n, , drop = FALSE]) %*%
                                          x[1:(n-k), , drop = FALSE]) %*% gamma0 ) )
  }
  m <- 1
  if( all(rhoK <= cn*sqrt(log(n)/n)) ){
    return(2*m)
  }else{
    for (k in (Kn+2):(n-1)) {
      m <- m+1
      tmpK <- rhoK
      rhoK[1:(Kn-1)] <- tmpK[2:Kn]
      rhoK[Kn] <- mean( abs( gamma0 %*% (t(x[(k+1):n, , drop = FALSE]) %*%
                                           x[1:(n-k), , drop = FALSE]) %*% gamma0 ) )
      if( all(rhoK <= cn*sqrt(log(n)/n)) ){
        return(2*m)
        break
      }
    }  # for
    
  }
  
}







