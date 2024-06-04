#' @name SpecTest
#' @title Statistical inference for high-dimensional spectral density matrix
#' @description \code{SpecTest()} implements a new global test proposed in
#'  Chang, Jiang, McElroy and Shao (2023) for the following hypothesis testing problem: 
#' \deqn{H_0:f_{i,j}(\omega)=0 \mathrm{\ for\ any\ }(i,j)\in\mathcal{I}\mathrm{\ and\ }
#' \omega \in \mathcal{J}\mathrm{\ \ versus\ \ H_1:H_0\ is\ not\ true.} }
#' 
#' 
#' @param X \eqn{{\bf X} = \{{\bf x_1}, \dots , {\bf x}_n \}}, a \eqn{n\times
#'   p} sample matrix, where \eqn{n} is the sample size and \eqn{p} is the 
#'   dimension of \eqn{{\bf x}_t}.
#' @param B Bootstrap times for generating multivariate normal distributed 
#' random vectors in calculating the critical value. Default is \code{B} \eqn{=2000}.
#' @param flag_c Bandwidth \eqn{c} of the flat-top kernel for estimating 
#' \eqn{f_{i,j}(\omega)}, where \eqn{c\in(0,1]}. Default is \code{flag_c} \eqn{=0.8}.
#' @param J.set Set \eqn{\mathcal{J}} for frequencies, a vector, used to calculate the test statistic.
#' @param cross.indices Set \eqn{\mathcal{I}} for \eqn{(i,j)}, a matrix with 2 columns,
#' used to calculate the test statistic.
#' 
#' @return An object of class "hdtstest" is a list containing the following
#'   components:
#'
#'   \item{Stat}{Numerical value which represents the value of test statistic.}
#'   \item{pval}{Numerical value which represents the p-value of the test.}
#'   \item{cri95}{Numerical value which represents the critical value of the test
#'   at the significance level 0.05.}
#'   \item{method}{A character string indicating what method was performed.}
#' @references Chang, J., Jiang, Q., McElroy, T. & Shao, X. (2023). 
#' \emph{Statistical inference for high-dimensional spectral density matrix}.
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
#' Stat <- res$Stat
#' Pvalue <- res$p.value
#' CriVal <- res$cri95
#' @export
SpecTest <- function(X, J.set, cross.indices, B = 1000, flag_c = 0.8)
{
  p <- ncol(X)
  n <- nrow(X)
  X <- t(X)
  K <- length(J.set)
  r <- dim(cross.indices)[1]
  l.band <- max(round(log(r)/10),1)
  
  ## ------------------------------------
  ##  Part II: compute spectral estimates
  ## ------------------------------------
  GhatC <- CmpGammaC(X)  # List(n) : p*p
  Shat <- matrix(l.band, nrow=p, ncol=p)
  
  ## compute spectral density estimates
  x.spcC <- SpecEstC(GhatC,n,p,r,K,cross.indices,J.set,l.band,flag_c) # List(K):p*p
  ## p.d. modification from Politis
  epsilon <- n^(-1.5)
  for(k in 1:length(J.set))
  {
    decomp <- getGCDc(x.spcC[[k]],p) 
    Lower.mat <- decomp[[1]]
    Diag.mat <- Re(decomp[[2]])
    Diag.mat <- pmax(Diag.mat,epsilon)
    x.spcC[[k]] <- Lower.mat %*% diag(Diag.mat) %*% Conj(t(Lower.mat))
  }
  x.spc <- array(unlist(x.spcC), c(p,p,K)) 
  
  ## ------------------------------------------
  ##  Part III: computing T1 and T2 statistics
  ## ------------------------------------------
  
  spc.stack <- matrix(NA, nrow = r, ncol = K)
  D.mat <- matrix(NA, nrow = r, ncol = K)
  for(j in 1:r)
  {
    k1 <- cross.indices[j,1]
    k2 <- cross.indices[j,2]
    spc.stack[j,] <- x.spc[k1,k2,] / sqrt(Shat[k1,k2])
  }
  
  delta.non <- sqrt(n)*spc.stack
  T2 <- max((Mod(delta.non))^2)
  
  
  ## ----------------------------------------------
  ## Part IV: computing T1* and T2* distributions
  ## ----------------------------------------------
  
  n.tilde <- n - 2*l.band -1
  C.hatC <- CEst2C(X,GhatC,n.tilde,n,p,r,cross.indices,l.band) # List(n.tilde): r*(2ln+1)
  C.hat <- array(unlist(C.hatC), c(r,(2*l.band+1), n.tilde))
  
  # ===================================================================
  #  critical value for studentized test statistic based on C.hat
  # ===================================================================
  bn <- BandEstC(matrix(unlist(C.hatC), c(r*(2*l.band+1), n.tilde)), 
                 n.tilde, r, l.band, type=1)
  eta <- etaC(n,p,B,n.tilde,bn,type=1) # n.tilde*B
  
  xi <- matrix(NA, nrow = 2*K*r, ncol = B)    # 2Kr*B
  mag <- matrix(NA, nrow = K*r, ncol = B)
  
  for(j in 1:r)
  {
    # cat('j=',j,'\r')
    Acos.mat <- matrix(NA, nrow = K, ncol = 2*l.band+1)  # K*(2ln+1)
    Asin.mat <- matrix(NA, nrow = K, ncol = 2*l.band+1)
    for (h in -l.band:l.band) {
      Acos.mat[,h+l.band+1] <- TaperFlatC(h/l.band,flag_c)*cos(J.set*h)/sqrt(l.band)
      Asin.mat[,h+l.band+1] <- -1*TaperFlatC(h/l.band,flag_c)*sin(J.set*h)/sqrt(l.band)
    }
    
    xi.cos <- Acos.mat %*% C.hat[j,,] %*% eta/sqrt(n.tilde) # K*B
    xi.sin <- Asin.mat %*% C.hat[j,,] %*% eta/sqrt(n.tilde)
    xi[((2*j-2)*K+1):((2*j-1)*K),] <- xi.cos
    xi[((2*j-1)*K+1):(2*j*K),] <- xi.sin
    mag[((j-1)*K+1):(j*K),] <- xi.cos^2 + xi.sin^2
  }
  
  
  T2.stars <- colMax(mag)
  names(T2) <- "Statistic"
  cri95=sort(T2.stars)[floor(.95*B)]
  names(cri95) <- "the critical value of the test at the significance level 0.05"
  METHOD = "Statistical inference for high-dimensional spectral density matrix"
  temp <- structure(list(statistic=T2, p.value=sum(T2.stars>T2)/B, 
               cri95=cri95,
               method = METHOD),
               class = "hdtstest")
  return(temp)
} 



#' @name SpecMulTest
#' @title Statistical inference for high-dimensional spectral density matrix
#' @description \code{SpecMulTest()} implements a new multiple test proposed in
#'  Chang, Jiang, McElroy and Shao (2023) for the \eqn{Q} hypothesis testing problems: 
#' \deqn{H_{0,q}:f_{i,j}(\omega)=0\mathrm{\ for\ any\ }(i,j)\in\mathcal{I}^{(q)}\mathrm{\ and\ }
#' \omega\in\mathcal{J}^{(q)}\mathrm{\ \ versus\ \ }
#' H_{1,q}:H_{0,q}\mathrm{\ is\ not\ true.} }
#' for \eqn{q\in\{1,\dots,Q\}}. 
#' 
#' 
#' @param Q Number of the hypothesis tests. 
#' @param PVal P-values for the \eqn{Q} hypothesis tests, a \eqn{Q} vector.
#' @param alpha The prescribed significance level. Default is 0.05.
#' @param seq_len Length used to take discrete points between 0 and 
#' \eqn{\sqrt(2\times\log(Q)-2\times\log(\log(Q))}. Default is 0.01.
#' 
#' 
#' @return An object of class "hdtstest" is a list containing the following
#'   components:
#'   \item{MultiTest}{Logical vector with length Q. If the element is \code{TRUE}, 
#'   it means rejecting the corresponding sub-null hypothesis, otherwise it means 
#'   not rejecting the corresponding sub-null hypothesis.}
#'   
#' @references Chang, J., Jiang, Q., McElroy, T. & Shao, X. (2023). 
#' \emph{Statistical inference for high-dimensional spectral density matrix}.
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
  
  #eFDR <- sum(MultiTest[,-H1_ind])/max(sum(MultiTest),1) 
  #epower <- sum(MultiTest[,H1_ind])/length(H1_ind)
  
  METHOD = "Statistical inference for high-dimensional spectral density matrix"
  
  temp <- structure(list(MultiTest=MultiTest,
                         method = METHOD),
                    class = "hdtstest")
  
  return(temp)
}



getGCDc <- function(Sigma,Rank)
{
  
  #############################
  #	getGCDc
  #		by Tucker McElroy
  #	Gets Generalized Cholesky Decomposition of complex Sigma,
  #	 allowing for zero Schur complements
  #      Rank is the presumed rank of the matrix, less than or equal
  #	  to the dimension of the matrix
  #	Output consists of the lower cholesky matrix,
  #	  and the diagonal matrix, of reduced dimension
  #
  #############################
  
  N <- dim(Sigma)[1]
  L.mat <- matrix(1,1,1)
  L.mat.inv <- L.mat
  D.mat <- Sigma[1,1]
  if(N > 1) {
    for(j in 2:N)
    {
      D.inv <- 1/D.mat
      D.inv[D.mat==0] <- 0
      new.sigma <- Sigma[j,1:(j-1)]
      if(j==2) { new.sigma <- as.matrix(new.sigma); L.mat <- as.matrix(L.mat) }
      new.l <- new.sigma %*% t(L.mat.inv)*D.inv
      new.l.tilde <- new.l %*% L.mat.inv
      L.mat <- cbind(L.mat,rep(0,(j-1)))
      L.mat <- rbind(L.mat,c(new.l,1))
      L.mat.inv <- cbind(L.mat.inv,rep(0,j-1))
      L.mat.inv <- rbind(L.mat.inv,c(-1*new.l.tilde,1))
      if(j==2) new.d <- Sigma[2,2] - Mod(new.l)^2*D.mat
      if(j > 2) new.d <- Sigma[j,j] - new.l %*% diag(D.mat) %*% t(Conj(new.l))
      new.d <- Re(new.d)
      if(new.d <= 0) { new.d <- 0 }
      D.mat <- c(D.mat,new.d)
    } }
  
  rank.index <- rank(D.mat,ties.method="first")
  dims <- seq(1,N)[rank.index > (N-Rank)]
  
  L.mat <- matrix(L.mat[,dims],nrow=N,ncol=length(dims))
  D.mat <- D.mat[dims]
  
  #	print(Lmat)
  #	print(Dmat)
  #	print(Lmat %*% diag(Dmat,nrow=length(dims)) %*% t(Lmat))
  
  return(list(L.mat,D.mat))
  #	return(list(L.mat,D.mat,dims))
}


getGCD <- function(Sigma,Rank)
{
  
  #############################
  #	getGCD
  #		by Tucker McElroy
  #	Gets Generalized Cholesky Decomposition of Sigma,
  #	 allowing for zero Schur complements
  #      Rank is the presumed rank of the matrix, less than or equal
  #	  to the dimension of the matrix
  #	Output consists of the lower cholesky matrix,
  #	  and the diagonal matrix, of reduced dimension
  #
  #############################
  
  N <- dim(Sigma)[1]
  L.mat <- matrix(1,1,1)
  L.mat.inv <- L.mat
  D.mat <- Sigma[1,1]
  if(N > 1) {
    for(j in 2:N)
    {
      
      D.inv <- 1/D.mat
      D.inv[D.mat==0] <- 0
      new.sigma <- Sigma[j,1:(j-1)]
      if(j==2) { new.sigma <- as.matrix(new.sigma); L.mat <- as.matrix(L.mat) }
      new.l <- new.sigma %*% t(L.mat.inv)*D.inv
      new.l.tilde <- new.l %*% L.mat.inv
      L.mat <- cbind(L.mat,rep(0,(j-1)))
      L.mat <- rbind(L.mat,c(new.l,1))
      L.mat.inv <- cbind(L.mat.inv,rep(0,j-1))
      L.mat.inv <- rbind(L.mat.inv,c(-1*new.l.tilde,1))
      if(j==2) new.d <- Sigma[2,2] - new.l^2*D.mat
      if(j > 2) new.d <- Sigma[j,j] - new.l %*% diag(D.mat) %*% t(new.l)
      if(new.d <= 0) { new.d <- 0 }
      D.mat <- c(D.mat,new.d)
    } }
  
  rank.index <- rank(D.mat,ties.method="first")
  dims <- seq(1,N)[rank.index > (N-Rank)]
  
  L.mat <- matrix(L.mat[,dims],nrow=N,ncol=length(dims))
  D.mat <- D.mat[dims]
  
  
  return(list(L.mat,D.mat))
  #	return(list(L.mat,D.mat,dims))
}



# index map
chi.map <- function(index,p) 
{ 
  ind.mat <- 0*diag(p)
  ind.mat[lower.tri(ind.mat)] <- seq(1,choose(p,2))
  return(which(ind.mat==index,arr.ind=TRUE))
}
chi.invmap <- function(indices,p)
{
  ind.mat <- 0*diag(p)
  ind.mat[lower.tri(ind.mat)] <- seq(1,choose(p,2))
  return(ind.mat[indices[1],indices[2]])
}

colMax <- function(data) apply(data, 2, max)







