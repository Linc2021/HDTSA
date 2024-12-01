#' @name SpecTest
#' @title Global testing for spectral density matrix
#' @description 
#' \code{SpecTest()} implements a new global test proposed in 
#' Chang et al. (2022) for the following hypothesis testing problem: 
#' \deqn{H_0:f_{i,j}(\omega)=0 \mathrm{\ for\ any\ }(i,j)\in \mathcal{I}\mathrm{\ and\ }
#' \omega \in \mathcal{J}\mathrm{\ \ versus\ \ }H_1:H_0\mathrm{\ is\ not\ true }\,,}
#' where \eqn{f_{i,j}(\omega)} represents the cross-spectral density between 
#' \eqn{ x_{t,i}} and \eqn{ x_{t,j}} at frequency \eqn{\omega} with \eqn{x_{t,i}} being 
#' the \eqn{i}-th component of the \eqn{p \times 1} times series \eqn{{\bf x}_t}.
#' Here, \eqn{\mathcal{I}} is the set of index pairs of interest, and
#' \eqn{\mathcal{J}} is the set of frequencies. 
#' 
#' 
#' 
#' @param X An \eqn{n\times p} data matrix \eqn{{\bf X} = ({\bf x}_1, \dots , {\bf x}_n)'},
#' where \eqn{n} is the number of observations of the \eqn{p\times 1} time
#' series \eqn{\{{\bf x}_t\}_{t=1}^n}.
#' @param B The number of bootstrap replications for generating multivariate normally
#' distributed random vectors when calculating the critical value. The default is 2000.
#' @param flag_c The bandwidth \eqn{c\in(0,1]} of the flat-top kernel for estimating 
#' \eqn{f_{i,j}(\omega)} [See (2) in Chang et al. (2022)]. The default is 0.8. 
#' @param J.set A vector representing the set \eqn{\mathcal{J}}
#' of frequencies.
#' @param cross.indices An \eqn{r \times 2} matrix representing the set
#' \eqn{\mathcal{I}} of \eqn{r} index pairs, where each row represents an index pair.
#' @seealso \code{\link{SpecMulTest}}
#' 
#' @return An object of class \code{"hdtstest"}, which contains the following
#'   components:
#'
#'   \item{Stat}{The test statistic of the test.}
#'   \item{pval}{The p-value of the test.}
#'   \item{cri95}{The critical value of the test
#'   at the significance level 0.05.}
#' @references Chang, J., Jiang, Q., McElroy, T. S., & Shao, X. (2022). 
#' Statistical inference for high-dimensional spectral density matrix.
#' \emph{arXiv preprint}. \doi{doi:10.48550/arXiv.2212.13686}.
#' 
#' @examples
#' # Example 1
#' ## Generate xt
#' n <- 200
#' p <- 10
#' flag_c <- 0.8
#' B <- 1000
#' burn <- 1000
#' z.sim <- matrix(rnorm((n+burn)*p),p,n+burn)
#' phi.mat <- 0.4*diag(p)
#' x.sim <- phi.mat %*% z.sim[,(burn+1):(burn+n)]
#' x <- x.sim - rowMeans(x.sim)
#' 
#' ## Generate the sets I and J
#' cross.indices <- matrix(c(1,2), ncol=2)
#' J.set <- 2*pi*seq(0,3)/4 - pi
#' 
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
    # METHOD = "Statistical inference for high-dimensional spectral density matrix"
    structure(list(statistic=T2, p.value=sum(T2.stars>T2)/B, 
                 cri95=cri95),
                 class = "hdtstest")
  }else{
    warning("n<=2ln+1 which does not satisfy the requirement of this test.")
    structure(list(statistic = NA, p.value = NA, 
                   cri95 = NA),
              class = "hdtstest")
  }
}



#' @name SpecMulTest
#' @title Multiple testing with FDR control for spectral density matrix
#' @description 
#' \code{SpecMulTest()} implements a new multiple testing procedure proposed in
#'  Chang et al. (2022) for the following \eqn{Q} hypothesis testing problems: 
#' \deqn{H_{0,q}:f_{i,j}(\omega)=0\mathrm{\ for\ any\ }(i,j)\in\mathcal{I}^{(q)}\mathrm{\ and\ }
#' \omega\in\mathcal{J}^{(q)}\mathrm{\ \ versus\ \ }
#' H_{1,q}:H_{0,q}\mathrm{\ is\ not\ true} }
#' for \eqn{q=1,\dots,Q}. 
#' Here, \eqn{f_{i,j}(\omega)} represents the cross-spectral density between 
#' \eqn{ x_{t,i}} and \eqn{ x_{t,j}} at frequency \eqn{\omega} with \eqn{x_{t,i}} being 
#' the \eqn{i}-th component of the \eqn{p \times 1} times series \eqn{{\bf x}_t},
#' and \eqn{\mathcal{I}^{(q)}} and \eqn{\mathcal{J}^{(q)}} refer to
#' the set of index pairs and the set of frequencies associated with the
#' \eqn{q}-th test, respectively.
#' 
#' 
#' @param Q The number of the hypothesis tests. 
#' @param PVal A vector of length \eqn{Q} representing p-values of the \eqn{Q} hypothesis tests.
#' @param alpha The prescribed level for the FDR control. The default is 0.05.
#' @param seq_len The step size for generating a sequence from 0 to
#' \eqn{\sqrt{2\times\log Q-2\times\log(\log Q )}}. The default is 0.01.
#' @seealso \code{\link{SpecTest}}
#' 
#' @return An object of class \code{"hdtstest"}, which contains the following
#'   component:
#'   \item{MultiTest}{A logical vector of length \eqn{Q}. If its \eqn{q}-th element is \code{TRUE}, 
#'   it indicates that \eqn{H_{0,q}} should be rejected. Otherwise,
#'   \eqn{H_{0,q}} should not be rejected.}
#'   
#' @references Chang, J., Jiang, Q., McElroy, T. S., & Shao, X. (2022). 
#' Statistical inference for high-dimensional spectral density matrix.
#' \emph{arXiv preprint}. \doi{doi:10.48550/arXiv.2212.13686}.
#' @examples
#' # Example 1
#' ## Generate xt
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
#' 
#' ## Generate the sets Iq and Jq
#' ISET <- list()
#' ISET[[1]] <- matrix(c(1,2),ncol=2)
#' ISET[[2]] <- matrix(c(1,3),ncol=2)
#' ISET[[3]] <- matrix(c(1,4),ncol=2)
#' ISET[[4]] <- matrix(c(1,5),ncol=2)
#' JSET <- as.list(2*pi*seq(0,3)/4 - pi)
#' 
#' ## Calculate Q p-values
#' PVal <- rep(NA,Q)
#' for (q in 1:Q) {
#'   cross.indices <- ISET[[q]]
#'   J.set <- JSET[[q]]
#'   temp.q <- SpecTest(t(x), J.set, cross.indices, B, flag_c)
#'   PVal[q] <- temp.q$p.value
#' }
#' res <- SpecMulTest(Q, PVal)
#' res
#' 
#' @useDynLib HDTSA
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @import Rcpp
#' @export
SpecMulTest <- function(Q, PVal, alpha=0.05, seq_len=0.01){
  
  Tseq = seq(0,sqrt(2*log(Q)-2*log(log(Q))), seq_len)
  
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
  
  # METHOD = "Statistical inference for high-dimensional spectral density matrix"
  
  structure(list(MultiTest=MultiTest),
                    class = "hdtstest")
}



# selection of ln
ChooseLn <- function(x, cn, Kn){
  n <- nrow(x)
  p <- ncol(x)
  gamma0 <- diag((diag(MatMult(t(x), x)))^(-1/2),nrow=p)
  rhoK <- rep(NA, Kn)
  for (k in 2:(Kn+1)) {
    rhoK[k-1] <- mean( abs( MatMult(MatMult(gamma0, MatMult(t(x[(k+1):n, , drop = FALSE]),
                                          x[1:(n-k), , drop = FALSE])), gamma0) ) )
  }
  m <- 1
  if( all(rhoK <= cn*sqrt(log(n)/n)) ){
    return(2*m)
  }else{
    for (k in (Kn+2):(n-1)) {
      m <- m+1
      tmpK <- rhoK
      rhoK[1:(Kn-1)] <- tmpK[2:Kn]
      rhoK[Kn] <- mean( abs( MatMult(MatMult(gamma0, MatMult(t(x[(k+1):n, , drop = FALSE]),
                                           x[1:(n-k), , drop = FALSE])), gamma0) ) )
      if( all(rhoK <= cn*sqrt(log(n)/n)) ){
        return(2*m)
        break
      }
    }  # for
    
  }
  
}







