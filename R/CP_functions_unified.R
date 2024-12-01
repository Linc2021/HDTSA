#' @title Estimating the matrix time series CP-factor model
#' @description \code{CP_MTS()} deals with the estimation of the CP-factor model for matrix time series:
#' \deqn{{\bf{Y}}_t = {\bf A \bf X}_t{\bf B}' +
#' {\boldsymbol{\epsilon}}_t, } where \eqn{{\bf X}_t = {\rm diag}(x_{t,1},\ldots,x_{t,d})} is a \eqn{d \times d}
#' unobservable diagonal matrix, \eqn{ {\boldsymbol{\epsilon}}_t }
#'  is a \eqn{p \times q} matrix white noise, \eqn{{\bf A}} and \eqn{{\bf B}} are, respectively, \eqn{p
#' \times d} and \eqn{q \times d} unknown constant matrices with their columns being
#' unit vectors, and \eqn{1\leq d < \min(p,q)} is an unknown integer.
#' Let \eqn{{\rm rank}(\mathbf{A}) = d_1}
#' and \eqn{{\rm rank}(\mathbf{B}) = d_2} with some unknown \eqn{d_1,d_2\leq d}.
#' This function aims to estimate \eqn{d, d_1, d_2} and the loading
#' matrices \eqn{{\bf A}} and \eqn{{\bf B}} using the methods proposed in Chang
#' et al. (2023) and Chang et al. (2024).
#' 
#' @details 
#' All three CP-decomposition methods involve the estimation of the autocovariance of
#' \eqn{ {\bf Y}_t} and \eqn{\xi_t} at lag \eqn{k}, which is defined as follows:
#' \deqn{\hat{\bf \Sigma}_{k} = T_{\delta_1}\{\hat{\boldsymbol{\Sigma}}_{\mathbf{Y},
#'  \xi}(k)\}\ \ {\rm with}\ \ \hat{\boldsymbol{\Sigma}}_{\mathbf{Y}, \xi}(k) = \frac{1}{n-k}
#' \sum_{t=k+1}^n(\mathbf{Y}_t-\bar{\mathbf{Y}})(\xi_{t-k}-\bar{\xi})\,,}
#' where \eqn{\bar{\bf Y} = n^{-1}\sum_{t=1}^n {\bf Y}_t}, \eqn{\bar{\xi}=n^{-1}\sum_{t=1}^n \xi_t}
#' and \eqn{T_{\delta_1}(\cdot)} is a threshold operator defined as
#' \eqn{T_{\delta_1}({\bf W}) = \{w_{i,j}1(|w_{i,j}|\geq \delta_1)\}} for any matrix
#' \eqn{{\bf W}=(w_{i,j})}, with the threshold level \eqn{\delta_1 \geq 0} and \eqn{1(\cdot)}
#' representing the indicator function. Chang et al. (2023) and Chang et al. (2024) suggest to choose
#' \eqn{\delta_1 = 0} when \eqn{p, q} are fixed and \eqn{\delta_1>0} when \eqn{pq \gg n}.
#' 
#' The refined estimation method involves
#' \deqn{\check{\bf \Sigma}_{k} =
#' T_{\delta_2}\{\hat{\mathbf{\Sigma}}_{\check{\mathbf{Y}}}(k)\}\ \ {\rm with}
#' \ \ \hat{\mathbf{\Sigma}}_{\check{\mathbf{Y}}}(k)=\frac{1}{n-k}
#' \sum_{t=k+1}^n(\mathbf{Y}_t-\bar{\mathbf{Y}}) \otimes {\rm vec}
#' (\mathbf{Y}_{t-k}-\bar{\mathbf{Y}})\,,}
#' where \eqn{T_{\delta_2}(\cdot)} is a threshold operator with the threshold level
#' \eqn{\delta_2 \geq 0}, and \eqn{{\rm vec}(\cdot)} is a vecterization operator
#' with \eqn{{\rm vec}({\bf H})} being the \eqn{(m_1m_2)\times 1} vector obtained by stacking
#' the columns of the \eqn{m_1 \times m_2} matrix \eqn{{\bf H}}. See Section 3.2.2 of Chang
#' et al. (2023) for details.
#' 
#' The unified estimation method involves
#' \deqn{\vec{\bf \Sigma}_{k}=
#' T_{\delta_3}\{\hat{\boldsymbol{\Sigma}}_{\vec{\mathbf{Y}}}(k)\}
#' \ \ {\rm with}\ \ \hat{\boldsymbol{\Sigma}}_{\vec{\mathbf{Y}}}(k)=\frac{1}{n-k}
#' \sum_{t=k+1}^n{\rm vec}({\mathbf{Y}}_t-\bar{\mathbf{Y}})\{{\rm vec}
#' (\mathbf{Y}_{t-k}-\bar{\mathbf{Y}})\}'\,,}
#' where \eqn{T_{\delta_3}(\cdot)} is a threshold operator with the threshold level
#' \eqn{\delta_3 \geq 0}. See Section 4.2 of Chang et al. (2024) for details.
#'
#'
#' @param Y An \eqn{n \times p \times q} array, where \eqn{n} is the number
#' of observations of the \eqn{p \times q} matrix time series \eqn{\{{\bf Y}_t\}_{t=1}^n}.
#' @param xi An \eqn{n \times 1} vector \eqn{\boldsymbol{\xi} = (\xi_1,\ldots, \xi_n)'},
#' where \eqn{\xi_t} represents a linear combination of \eqn{{\bf Y}_t}.
#' If \code{xi = NULL} (the default), \eqn{\xi_{t}} is determined by the PCA
#' method introduced in Section 5.1 of Chang et al. (2023). Otherwise, \code{xi}
#' can be given by the users.
#' @param Rank A list containing the following components: \code{d} representing
#' the number of columns of \eqn{{\bf A}} and \eqn{{\bf B}}, \code{d1} representing
#'  the rank of \eqn{{\bf A}}, and \code{d2} representing the rank of \eqn{{\bf B}}.
#' If set to \code{NULL} (default), \eqn{d}, \eqn{d_1}, and \eqn{d_2} will be estimated.
#'  Otherwise, they can be given by the users.
#' @param lag.k The time lag \eqn{K} used to calculate the nonnegative definite 
#' matrices \eqn{\hat{\mathbf{M}}_1} and \eqn{\hat{\mathbf{M}}_2} when \code{method = "CP.Refined"}
#'  or \code{method = "CP.Unified"}:
#'  \deqn{\hat{\mathbf{M}}_1\ =\
#'   \sum_{k=1}^{K} \hat{\bf \Sigma}_{k} \hat{\bf \Sigma}_{k}'\ \ {\rm and}
#'   \ \ \hat{\mathbf{M}}_2\ =\ \sum_{k=1}^{K} \hat{\bf \Sigma}_{k}' \hat{\bf \Sigma}_{k}\,,
#'   }
#'   where \eqn{\hat{\bf \Sigma}_{k}} is an estimate of the cross-covariance between
#'   \eqn{ {\bf Y}_t} and \eqn{\xi_t} at lag \eqn{k}. See 'Details'. The default is 20.
#' @param lag.ktilde The time lag \eqn{\tilde K} involved in the unified
#' estimation method [See (16) in Chang et al. (2024)], which is used
#' when \code{method = "CP.Unified"}. The default is 10.
#' @param method A string indicating which CP-decomposition method is used. Available options include:
#'  \code{"CP.Direct"} (the default) for the direct estimation method
#'  [See Section 3.1 of Chang et al. (2023)], \code{"CP.Refined"} for the refined estimation
#'  method [See Section 3.2 of Chang et al. (2023)], and \code{"CP.Unified"} for the
#'  unified estimation method [See Section 4 of Chang et al. (2024)].
#'  The validity of methods \code{"CP.Direct"} and \code{"CP.Refined"} depends on the assumption
#'  \eqn{d_1=d_2=d}. When \eqn{d_1,d_2 \leq d}, the method \code{"CP.Unified"} can be applied.
#'  See Chang et al. (2024) for details.
#'  
#' @param thresh1  Logical. If \code{FALSE} (the default), no thresholding will
#'   be applied in \eqn{\hat{\bf \Sigma}_{k}}, which indicates that the threshold level
#'  \eqn{\delta_1=0}. If \code{TRUE}, \eqn{\delta_1} will be set through \code{delta1}.
#'   \code{thresh1} is used for all three methods. See 'Details'.
#' @param thresh2  Logical. If \code{FALSE} (the default), no thresholding will
#'   be applied in \eqn{\check{\bf \Sigma}_{k}}, which indicates that the threshold level
#'   \eqn{\delta_2=0}. If \code{TRUE}, \eqn{\delta_2} will be set through \code{delta2}.
#'   \code{thresh2} is used only when \code{method = "CP.Refined"}. See 'Details'.
#' @param thresh3  Logical. If \code{FALSE} (the default), no thresholding will
#'   be applied in \eqn{\vec{\bf \Sigma}_{k}}, which indicates that the threshold level
#'   \eqn{\delta_3=0}. If \code{TRUE}, \eqn{\delta_3} will be set through \code{delta3}.
#'   \code{thresh3} is used only when \code{method = "CP.Unified"}. See 'Details'.
#' @param delta1  The value of the threshold level \eqn{\delta_1}. The default is
#'  \eqn{ \delta_1 = 2 \sqrt{n^{-1}\log (pq)}}.
#' @param delta2  The value of the threshold level \eqn{\delta_2}. The default is
#'  \eqn{ \delta_2 = 2 \sqrt{n^{-1}\log (pq)}}.
#' @param delta3  The value of the threshold level \eqn{\delta_3}. The default is
#'  \eqn{ \delta_3 = 2 \sqrt{n^{-1}\log(pq)}}.
#'   
#' @return An object of class \code{"mtscp"}, which contains the following
#'   components:
#'   \item{A}{The estimated \eqn{p \times \hat{d}} left loading matrix \eqn{\hat{\bf A}}.}
#'   \item{B}{The estimated \eqn{q \times \hat{d}} right loading matrix \eqn{\hat{\bf B}}.}
#'   \item{f}{The estimated latent process \eqn{\hat{x}_{t,1},\ldots,\hat{x}_{t,\hat{d}}}.}
#'   \item{Rank}{The estimated \eqn{\hat{d}_1,\hat{d}_2}, and \eqn{\hat{d}}.}
#'   \item{method}{A string indicating which CP-decomposition method is used.}
#'
#'
#' @references
#'   Chang, J., Du, Y., Huang, G., & Yao, Q. (2024). Identification and
#'  estimation for matrix time series CP-factor models. \emph{arXiv preprint}.
#'  \doi{doi:10.48550/arXiv.2410.05634}.
#'  
#'   Chang, J., He, J., Yang, L., & Yao, Q. (2023). Modelling matrix time series via a tensor CP-decomposition.
#'   \emph{Journal of the Royal Statistical Society Series B: Statistical Methodology}, \strong{85}, 127--148. 
#'   \doi{doi:10.1093/jrsssb/qkac011}.
#'   
#'
#' @examples
#' # Example 1.
#' p <- 10
#' q <- 10
#' n <- 400
#' d = d1 = d2 <- 3
#' ## DGP.CP() generates simulated data for the example in Chang et al. (2024).
#' data <- DGP.CP(n, p, q, d, d1, d2)
#' Y <- data$Y
#' 
#' ## d is unknown
#' res1 <- CP_MTS(Y, method = "CP.Direct")
#' res2 <- CP_MTS(Y, method = "CP.Refined")
#' res3 <- CP_MTS(Y, method = "CP.Unified")
#' 
#' ## d is known
#' res4 <- CP_MTS(Y, Rank = list(d = 3), method = "CP.Direct")
#' res5 <- CP_MTS(Y, Rank = list(d = 3), method = "CP.Refined")
#' 
#' 
#' # Example 2.
#' p <- 10
#' q <- 10
#' n <- 400
#' d1 = d2 <- 2
#' d <-3
#' data <- DGP.CP(n, p, q, d, d1, d2)
#' Y1 <- data$Y
#' 
#' ## d, d1 and d2 are unknown
#' res6 <- CP_MTS(Y1, method = "CP.Unified")
#' ## d, d1 and d2 are known
#' res7 <- CP_MTS(Y1, Rank = list(d = 3, d1 = 2, d2 = 2), method = "CP.Unified")
#' 
#' @export
#' @useDynLib HDTSA
#' @importFrom stats arima.sim rnorm runif

CP_MTS = function(Y, xi = NULL, Rank = NULL, lag.k = 20, lag.ktilde = 10,
                  method = c("CP.Direct","CP.Refined","CP.Unified"),
                  thresh1 = FALSE, thresh2 = FALSE, thresh3 = FALSE,
                  delta1 = 2 * sqrt(log(dim(Y)[2] * dim(Y)[3]) / dim(Y)[1]),
                  delta2 = delta1, delta3 = delta1){
  n <- dim(Y)[1]; p <- dim(Y)[2]; q <- dim(Y)[3];
  if(is.null(xi)){
    xi <- est.xi(Y)$xi
  }
  method <- match.arg(method)
  if(method == "CP.Direct"){
    S_yxi_1 <- Autocov_xi_Y(Y, xi, k = 1, thresh = thresh1, delta = delta1)
    S_yxi_2 <- Autocov_xi_Y(Y, xi, k = 2, thresh = thresh1, delta = delta1)
    if(p > q){
      ##(1) estimation of d
      K1 <- t(S_yxi_1) %*% S_yxi_1
      eg1 <- eigen(K1)
      w <- eg1$values
      ww <- w[-1] / w[-length(w)]
      d <- which(ww == min(ww[1:floor(0.75 * q)]))
      if(!is.null(Rank)){
        if(!is.null(Rank$d)){
          d <- Rank$d
        }
        else{stop("List Rank without d, use Rank=list(d=?)")}
      }
      if (d > 1){
        K1 <- eg1$vectors[, 1:d]%*%diag(eg1$values[1:d])%*%t(eg1$vectors[, 1:d]);
      }else{
        K1 <- eg1$vectors[, 1]%*%diag(eg1$values[1], 1)%*%t(eg1$vectors[, 1]);
      }
      K2 <- t(S_yxi_1) %*% S_yxi_2;
      
      ##(2) estimation of A and B
      Geg <- geigen::geigen(K2, K1);
      evalues <- Geg$values[which(Mod(Geg$values) <= 10^5 & Geg$values != 0)]
      
      Bl <- Geg$vectors[, which(Geg$values %in% evalues), drop = FALSE]
      A <- apply(S_yxi_1 %*% Bl, 2, l2s)
      Al <- t(MASS::ginv(A))
      B <- apply(t(S_yxi_1) %*% Al, 2, l2s)
    }else{
      
      ##(1) estimation of d
      K1 <- S_yxi_1 %*% t(S_yxi_1)
      eg1 <- eigen(K1)
      w <- eg1$values
      ww <- w[-1] / w[-length(w)]
      d <- which(ww == min(ww[1:floor(0.75 * p)]))
      if(!is.null(Rank)){
        if(!is.null(Rank$d)){
          d <- Rank$d
        }
        else{stop("List Rank without d, use Rank=list(d=?)")}
      }
      
      if (d > 1){
        K1 <- eg1$vectors[, 1:d] %*% diag(eg1$values[1:d]) %*% t(eg1$vectors[, 1:d]);
      }else{
        K1 <- eg1$vectors[, 1] %*% diag(eg1$values[1], 1) %*% t(eg1$vectors[, 1]);
      }
      K2 <- S_yxi_1 %*% t(S_yxi_2);
      ##(2) estimation of A and B
      Geg <- geigen::geigen(K2, K1);
      evalues <- Geg$values[which(Mod(Geg$values) <= 10^5 & Geg$values!=0)]
      Al <- Geg$vectors[, which(Geg$values %in% evalues), drop = FALSE]
      B <- apply(t(S_yxi_1) %*% Al, 2, l2s)
      Bl <- t(MASS::ginv(B))
      A <- apply(S_yxi_1 %*% Bl, 2, l2s)
    }
    ##(3) estimation of Xt
    H <- matrix(NA, p * q, d)
    for (ii in 1:d) {
      H[, ii] <- B[, ii] %x% A[, ii]
    }
    f <- Vec.tensor(Y) %*% H %*% MASS::ginv(t(H) %*% H)
    
    if(is.complex(A) == T || is.complex(B) == T ){
      A <- Complex2Real(A)
      B <- Complex2Real(B)
      f <- Complex2Real(f)
    }
    # METHOD <- c("Estimation of matrix CP-factor model",paste("Method:", method))
    con <- structure(list(A = A,B = B,f = f,Rank = list(d = d), method = method),
                     class = "mtscp")
    return(con)
  }
  if(method == "CP.Refined"){
    ##(1) estimation of P,Q and d
    M1 = M2 <- 0
    dmax <- round(min(p, q) * 0.75)
    
    for (kk in 1:lag.k){
      S_yxi_k <- Autocov_xi_Y(Y, xi, k = kk, thresh = thresh1, delta = delta1)
      M1 <- M1 + S_yxi_k %*% t(S_yxi_k)
      M2 <- M2 + t(S_yxi_k) %*% S_yxi_k
    }
    ev_M1 <- eigen(M1)
    ev_M2 <- eigen(M2)
    
    d1 <-  which.max(ev_M1$values[1:dmax] / ev_M1$values[2:(dmax + 1)])
    d2 <-  which.max(ev_M2$values[1:dmax] / ev_M2$values[2:(dmax + 1)])
    
    d <- ifelse(p > q, d1, d2)
    
    if(!is.null(Rank)){
      if(!is.null(Rank$d)){
        d <- Rank$d
      }
      else{stop("List Rank without d, use Rank=list(d=?)")}
    }
    P <- ev_M1$vectors[, 1:d, drop = FALSE]
    Q <- ev_M2$vectors[, 1:d, drop = FALSE]
    
    if(d == 1){
      A <- as.matrix(P)
      B <- as.matrix(Q)
      
      f <- vector()
      for (tt in 1:n) {
        f[tt] <- t(A) %*% Y[tt, , ] %*% B
      }
      f <- as.matrix(f)
      
    }else{
      ##(2) estimation of U and V
      Z <- array(NA, dim = c(n, d, d))
      for (tt in 1:n) {
        Z[tt, , ] <- t(P) %*% Y[tt, , ] %*% Q
      }
      xi <- est.xi(Z)
      if(thresh2){
        w_hat <- xi$w_hat
        Xi <- diag(1, p) %x% ((Q %x% P) %*% as.matrix(w_hat))
        sigma_ycheck_1 <- Sigma_Ycheck(Y, 1)
        sigma_ycheck_1 <- thresh_C(sigma_ycheck_1, delta2)
        sigma_ycheck_2 <- Sigma_Ycheck(Y, 2)
        sigma_ycheck_2 <- thresh_C(sigma_ycheck_2, delta2)
        S_Zxi_1 <- t(P) %*% t(Xi) %*% sigma_ycheck_1 %*% Q
        S_Zxi_2 <- t(P) %*% t(Xi) %*% sigma_ycheck_2 %*% Q
      }
      else{
        S_Zxi_1 <- Autocov_xi_Y(Z, xi$xi, k = 1)
        S_Zxi_2 <- Autocov_xi_Y(Z, xi$xi, k = 2)
      }
      
      vl <- eigen(MASS::ginv(t(S_Zxi_1) %*% S_Zxi_1) %*% t(S_Zxi_1) %*% S_Zxi_2)$vectors ##MASS
      ul <- eigen(MASS::ginv(S_Zxi_1 %*% t(S_Zxi_1)) %*% S_Zxi_1 %*% t(S_Zxi_2))$vectors
      
      U <- apply(S_Zxi_1 %*% vl, 2, l2s)
      V <- apply(t(S_Zxi_1) %*% ul, 2, l2s)
      
      ##(3) estimation of A and B
      A <- P %*% U
      B <- Q %*% V
      
      ##(4) estimation of Xt
      W <- matrix(NA, d^2, d)
      
      for (ii in 1:d) {
        W[, ii] <- V[, ii] %x% U[, ii]
      }
      
      f <- Vec.tensor(Z) %*% W %*% solve(t(W) %*% W)
      
      if(is.complex(A) == T || is.complex(B) == T ){
        A <- Complex2Real(A)
        B <- Complex2Real(B)
        f <- Complex2Real(f)
      }
    }
    # METHOD <- c("Estimation of matrix CP-factor model",paste("Method:", method))
    con <- structure(list(A = A,B = B,f = f,Rank = list(d = d), method = method),
                     class = "mtscp")
    return(con)
  }
  if(method == "CP.Unified"){
    ##(1) estimation of P,Q and d1,d2
    if(is.null(Rank)){
      
      PQ_hat_tol <- est.d1d2.PQ(Y, xi, K = lag.k, thresh = thresh1, delta = delta1)
      
      d1 <- PQ_hat_tol$d1_hat
      d2 <- PQ_hat_tol$d2_hat
      d  <- NULL
      P  <- PQ_hat_tol$P_hat
      Q  <- PQ_hat_tol$Q_hat
      
      if(d1 == 1 || d2 == 1){d <- d1 * d2}
      
    }else{
      if(all(!is.null(Rank$d1), !is.null(Rank$d1), !is.null(Rank$d2))){
        d  <- Rank$d
        d1 <- Rank$d1
        d2 <- Rank$d2
      }
      else{stop("List Rank without d, d1 and d2, use Rank=list(d=?, d1=?, d2=?)")}
      
      PQ_hat_tol <- est.PQ(Y, xi, d1, d2, K = lag.k, thresh = thresh1, delta = delta1)
      
      P   <- PQ_hat_tol$P_hat
      Q   <- PQ_hat_tol$Q_hat
      
    }
    ##(2) estimation of W* = (v1*u1,v2*u2,...,vd*ud)H = WH
    if(d1 == 1 & d2 == 1){
      d <- 1
      f <- vector()
      for (tt in 1:n) {
        f[tt] = t(P) %*% Y[tt, , ] %*% Q
      }
      A <- P
      B <- Q
      
      # METHOD <- c("Estimation of matrix CP-factor model",paste("Method:",method))
      rank <- list(d = d, d1 = d1, d2 = d2)
      con <- structure(list(A = as.matrix(A), B = as.matrix(B), 
                            f = as.matrix(f), Rank = rank, method = method),
                       class = "mtscp")
      
      return(con)
    }else{
      if(is.null(d)){
        W_hat_tol  <-  est.d.Wf(Y, P, Q, Ktilde = lag.ktilde, thresh = thresh3, delta = delta3)
        d          <-  W_hat_tol$d_hat
        W          <-  W_hat_tol$W_hat
        f          <-  W_hat_tol$f_hat
      }else{
        d          <-  d
        W_hat_tol  <-  est.Wf(Y, P, Q, d, Ktilde = lag.ktilde, thresh = thresh3, delta = delta3)
        W          <-  W_hat_tol$W_hat
        f          <-  W_hat_tol$f_hat
      }
      
      ##(3) estimation of U and V
      if(d1 == 1 || d2 == 1){
        Theta <- NULL
        if(d1 == 1){
          U <- 1;
          V <- W;
          stop("d1 equal to 1, V cannot be identified uniquely!")
        }
        if(d2 == 1){
          U <- W;
          V <- 1;
          stop("d2 equal to 1, U cannot be identified uniquely!")
        }
        if(d1 == 1 & d2 == 1){
          U <- 1;
          V <- 1;
        }
        U <- as.matrix(U)
        V <- as.matrix(V)
      }else{
        
        UV_hat_tol <- est.UV.JAD(W, d1, d2, d)
        
        U          <- UV_hat_tol$U
        V          <- UV_hat_tol$V
        Theta      <- UV_hat_tol$Theta
      }
      
      
      ##(4) estimation of A and B
      A <- P %*% U
      B <- Q %*% V
      # METHOD <- c("Estimation of matrix CP-factor model",paste("Method:",method))
      rank <- list(d = d, d1 = d1, d2 = d2)
      con <- structure(list(A = as.matrix(A), B = as.matrix(B),
                            f = as.matrix(f), Rank = rank, method = method),
                       class = "mtscp")
      
      return(con)
    }
    
  }
}

rho2.loss = function(A_hat,A){
  max(apply(1-(t(A_hat)%*%A)^2,2,min))
}

l2s = function(x){x/sqrt(sum(x^2))}

fnorm = function(x){sqrt(sum(x^2))}

Complex2Real = function(A){
  REA = round(Re(A),8)
  IMA = round(Im(A),8)
  
  real.index  =  which(IMA[1,] == 0)
  
  if(length(real.index) == 0){
    complex_real  =  REA
    complex_image =  IMA
    
    complex_take  = which(duplicated(complex_real[1,]) == T)
    
    real = as.matrix(complex_real[,complex_take])
    
    img  = as.matrix(complex_image[,complex_take])
    
    new_A = cbind(real,img)
  }else{
    real.vector   = REA[,real.index]
    
    complex_real  =  REA[,-real.index]
    complex_image =  IMA[,-real.index]
    
    complex_take  = which(duplicated(complex_real[1,]) == T)
    
    real = as.matrix(complex_real[,complex_take])
    
    img  = as.matrix(complex_image[,complex_take])
    
    new_A = cbind(real.vector,real,img)
    
    colnames(new_A) = NULL
  }
  
  return(new_A)
}


Vec.tensor = function(Y){
  p = dim(Y)[2];q = dim(Y)[3];
  
  if(p == q & q == 1){
    Y_tilde = apply(Y,1,FUN = as.vector)
  }else{
    Y_tilde = t(apply(Y,1,FUN = as.vector))
  }
  return(Y_tilde)
  
}


#' @title Generating simulated data for the example in Chang et al. (2024)
#' @description \code{DGP.CP()} function generates simulated data following the 
#' data generating process described in Section 7.1 of Chang et al. (2024). 
#'  
#'
#' @param n  Integer. The number of observations of the \eqn{p \times q} matrix 
#' time series \eqn{{\bf Y}_t}.
#' @param p  Integer. The number of rows of \eqn{{\bf Y}_t}.
#' @param q  Integer. The number of columns of \eqn{{\bf Y}_t}.
#' @param d  Integer. The number of columns of the factor loading matrices \eqn{\bf A} 
#' and \eqn{\bf B}.
#' @param d1  Integer. The rank of the \eqn{p \times d} matrix \eqn{\bf A}.
#' @param d2  Integer. The rank of the \eqn{q \times d} matrix \eqn{\bf B}.
#'
#' @seealso \code{\link{CP_MTS}}.
#' @return A list containing the following
#'   components:
#'   \item{Y}{An \eqn{n \times p \times q} array.}
#'   \item{A}{The \eqn{p \times d} left loading matrix \eqn{\bf A}.}
#'   \item{B}{The \eqn{q \times d} right loading matrix \eqn{\bf B}.}
#'   \item{X}{An \eqn{n \times d \times d} array.}
#' @references
#'   Chang, J., Du, Y., Huang, G., & Yao, Q. (2024). Identification and
#'  estimation for matrix time series CP-factor models. \emph{arXiv preprint}.
#'  \doi{doi:10.48550/arXiv.2410.05634}.
#'  
#' @details We generate
#' \deqn{{\bf{Y}}_t = {\bf A \bf X}_t{\bf B}' + {\boldsymbol{\epsilon}}_t }
#' for any \eqn{t=1, \ldots, n}, where \eqn{{\bf X}_t = {\rm diag}({\bf x}_t)} 
#' with \eqn{{\bf x}_t = (x_{t,1},\ldots,x_{t,d})'} being a \eqn{d \times 1} time series,
#' \eqn{ {\boldsymbol{\epsilon}}_t } is a \eqn{p \times q} matrix white noise,
#' and \eqn{{\bf A}} and \eqn{{\bf B}} are, respectively, \eqn{p\times d} and 
#' \eqn{q \times d} factor loading matrices. \eqn{\bf A}, \eqn{{\bf X}_t}, and \eqn{\bf B}
#' are generated based on the data generating process described in Section 7.1 of 
#' Chang et al. (2024) and satisfy \eqn{{\rm rank}({\bf A})=d_1} and 
#' \eqn{{\rm rank}({\bf B})=d_2}, \eqn{1 \le d_1, d_2 \le d}.
#'
#' @examples
#' p <- 10
#' q <- 10
#' n <- 400
#' d = d1 = d2 <- 3
#' data <- DGP.CP(n,p,q,d1,d2,d)
#' Y <- data$Y
#' 
#' ## The first observation: Y_1
#' Y[1, , ]
#' @export
DGP.CP = function(n,p,q,d,d1,d2){
  
  par_A = c(-3,3)
  par_B = c(-3,3)
  par_X = c(0.6,0.95)
  par_E = 1
  
  Input = list(n = n,p = p,q = q,d1 = d1,d2 = d2,d = d)
  
  A_inl = matrix(runif(p*d,par_A[1],par_A[2]),p,d)
  B_inl = matrix(runif(q*d,par_B[1],par_B[2]),q,d)
  
  svd_A = svd(A_inl)
  P = svd_A$u[,1:d1]
  
  A = (svd_A$u[,1:d1]) %*% diag(svd_A$d[1:d1],nrow = d1,ncol = d1) %*%t(svd_A$v[,1:d1])
  
  U = t(P)%*%apply(A,2,l2s)
  
  A_s = P%*%U
  
  
  svd_B = svd(B_inl)
  Q = svd_B$u[,1:d2]
  B = (svd_B$u[,1:d2]) %*% diag(svd_B$d[1:d2],nrow = d2,ncol = d2) %*% (t(svd_B$v[,1:d2]))
  V = t(Q)%*%apply(B,2,l2s)
  
  B_s = Q%*%V
  
  
  W = matrix(NA,d1*d2,d)
  
  for (ii in 1:d) {
    W[,ii] = V[,ii]%x%U[,ii]
  }
  
  X = array(0,c(n,d,d))
  X_m = matrix(NA,n,d)
  
  signal = runif(d,-1,1)
  signal[signal>=0] <-  1
  signal[signal< 0] <- -1
  par_ar = runif(d,par_X[1],par_X[2])
  for (ii in 1:d) {
    xd = arima.sim(model = list(ar = par_ar[ii]*signal[ii]),n = n)
    X_m[,ii] =  xd*(apply(A, 2, fnorm)*apply(B, 2, fnorm))[ii]
    X[,ii,ii] <- xd*(apply(A, 2, fnorm)*apply(B, 2, fnorm))[ii]
    
  }
  
  S_m = X_m%*%t(W)
  W_star = svd(S_m)$v[,1:d]
  
  Y = S = array(NA,dim = c(n,p,q))
  for (tt in 1:n) {
    S[tt,,] <- A_s%*%X[tt,,]%*%t(B_s)
    Y[tt,,] <- S[tt,,] + matrix(rnorm(p*q,0,par_E),p,q)
  }
  
  return(list(Y      =  Y,
              A      =  A_s,
              B      =  B_s,
              X      =  X
  ))
  
}


Autocov_xi_Y = function(Y, xi, k, thresh = FALSE, delta = NULL){
  
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  q <- dim(Y)[3]
  # k <- lag.k
  
  Y_mean <- 0
  xi_mean <- 0
  for (ii in 1:n) {
    Y_mean <- Y_mean + Y[ii,,]
    xi_mean <- xi_mean + xi[ii]
  }
  Y_mean  <- Y_mean/n
  xi_mean <- xi_mean/n
  
  Sigma_Y_xi_k <- 0
  for (ii in (k+1):n) {
    Sigma_Y_xi_k <- Sigma_Y_xi_k + (Y[ii,,] - Y_mean)*(xi[ii-k] - xi_mean)
  }
  Sigma_Y_xi_k <- Sigma_Y_xi_k/(n-k)
  if(thresh){
    Sigma_Y_xi_k <- thresh_C(Sigma_Y_xi_k, delta)
  }
  return(Sigma_Y_xi_k)
}

est.xi  = function(Y, thresh_per = 0.99, d_max = 20){
  
  n = dim(Y)[1];p = dim(Y)[2];q = dim(Y)[3];
  
  xi.mat = Vec.tensor(Y)
  
  if(n > p*q){
    eig_xi.mat = eigen(MatMult(t(xi.mat),xi.mat))
    cfr =  cumsum(eig_xi.mat$values)/sum(eig_xi.mat$values)
    d_hat = min(which(cfr > thresh_per))
    d_fin = min(d_max,d_hat)
    w_hat = eig_xi.mat$vectors[,1:d_fin, drop=FALSE]
    if(d_fin == 1){
      sign_value = adjust_sign(w_hat[, 1])
      w_hat = w_hat * sign_value
    }else{
      column_signs = apply(w_hat, 2, adjust_sign)
      w_hat = w_hat %*% diag(column_signs)
    }
    xi.f = xi.mat%*%w_hat
    xi   = rowMeans(xi.f)
    w_hat = rowMeans(w_hat)
  }else{
    eig_xi.mat = eigen(MatMult(xi.mat,t(xi.mat)))
    cfr = cumsum(eig_xi.mat$values)/sum(eig_xi.mat$values)
    d_hat = min(which(cfr > thresh_per))
    d_fin = min(d_max,d_hat)
    xi.f1 = as.matrix(eig_xi.mat$vectors[,1:d_fin, drop=FALSE])
    if(d_fin == 1){
      sign_value = adjust_sign(xi.f1[, 1])
      xi.f1 = xi.f1 * sign_value
    }else{
      column_signs = apply(xi.f1, 2, adjust_sign)
      xi.f1 = xi.f1 %*% diag(column_signs)
    }
    weight = sqrt(eig_xi.mat$values[1:d_fin])
    xi.f1 = xi.f1%*%diag(weight)
    xi = rowMeans(xi.f1)
    w_hat = as.matrix(t(t(xi.f1)%*%xi.mat))
    w_hat = rowMeans(w_hat)
  }
  return(list(xi=xi, w_hat = w_hat))
}

est.d1d2.PQ = function(Y,xi,K = 10, thresh = FALSE, delta = NULL){
  n = dim(Y)[1];p = dim(Y)[2];q = dim(Y)[3];
  
  d2_list = d1_list =vector()
  M1 = M2 = 0
  dmax = round(min(p,q)*0.75)
  P_list = Q_list = list()
  for (kk in 1:K){
    
    S_yxi_k = Autocov_xi_Y(Y,xi, k = kk, thresh = thresh, delta = delta)
    
    M1 = M1 + S_yxi_k%*%t(S_yxi_k)
    M2 = M2 + t(S_yxi_k)%*%S_yxi_k
    
    ev_M1 = eigen(M1)
    ev_M2 = eigen(M2)
    
    d1_list[kk] =  which.max(ev_M1$values[1:dmax]/ev_M1$values[2:(dmax+1)])
    d2_list[kk] =  which.max(ev_M2$values[1:dmax]/ev_M2$values[2:(dmax+1)])
    
    P_list[[kk]] = ev_M1$vectors[,1:(d1_list[kk]), drop = FALSE]
    Q_list[[kk]] = ev_M2$vectors[,1:(d2_list[kk]), drop = FALSE]
    
  }
  
  d1_list[1] = 0
  d2_list[1] = 0
  
  
  d1_hat =  d1_list[K]
  d2_hat =  d2_list[K]
  P_hat  =  P_list[[K]]
  Q_hat  =  Q_list[[K]]
  
  return(list(d1_hat  = d1_hat,
              d2_hat  = d2_hat,
              P_hat   = P_hat,
              Q_hat   = Q_hat,
              d1_list = d1_list,
              d2_list = d2_list,
              P_list  = P_list,
              Q_list  = Q_list))
}

est.PQ = function(Y,xi,d1,d2,K = 20, thresh = FALSE, delta = NULL){
  n = dim(Y)[1];p = dim(Y)[2];q = dim(Y)[3];
  
  M1 = M2 = 0
  dmax = round(min(p,q)*0.75)
  P_list = Q_list = list()
  for (kk in 1:K){
    
    S_yxi_k = Autocov_xi_Y(Y,xi, k = kk, thresh = thresh, delta = delta)
    
    M1 = M1 + S_yxi_k%*%t(S_yxi_k)
    M2 = M2 + t(S_yxi_k)%*%S_yxi_k
    
    ev_M1 = eigen(M1)
    ev_M2 = eigen(M2)
    
    P_list[[kk]] = ev_M1$vectors[ ,1:(d1), drop = FALSE]
    Q_list[[kk]] = ev_M2$vectors[ ,1:(d2), drop = FALSE]
    
  }
  
  P_hat  =  P_list[[K]]
  Q_hat  =  Q_list[[K]]
  
  return(list(P_hat   = P_hat,
              Q_hat   = Q_hat,
              P_list  = P_list,
              Q_list  = Q_list))
}


est.d.Wf = function(Y,P,Q, Ktilde = 10, thresh = FALSE, delta = NULL){
  n = dim(Y)[1];p = dim(Y)[2];q = dim(Y)[3];
  d1 = NCOL(P);d2 = NCOL(Q);
  
  Z = array(NA,dim = c(n,d1,d2))
  for (tt in 1:n) {
    Z[tt,,] = t(P)%*%Y[tt,,]%*%Q
  }
  
  Z_tilde = Vec.tensor(Z)
  
  M = 0
  W_list = f_list = list()
  d_list = vector()
  dmax = d1*d2
  dstar = max(d1,d2)
  
  for (kk in 1:Ktilde) {
    
    if(thresh){
      Y_2d <- Vec.tensor(Y)
      S_ytilde_k <- sigmak(t(Y_2d), as.matrix(colMeans(Y_2d)), n = n, k = kk)
      S_ytilde_k <- thresh_C(S_ytilde_k, delta)
      S_ztilde_k <- MatMult(MatMult(t(Q) %x% t(P), S_ytilde_k), Q %x% P)
    }
    else{
    S_ztilde_k = sigmak(t(Z_tilde), as.matrix(colMeans(Z_tilde)),n = n, k= kk)
    }
    
    M = M + S_ztilde_k%*%t(S_ztilde_k)
    
    ev_M = eigen(M)
    
    evalues = ev_M$values
    
    
    d_list[kk] =  max(which.max(evalues[1:(dmax-1)]/evalues[2:(dmax)]),dstar)
    
    W_list[[kk]] = ev_M$vectors[,1:(d_list[kk]), drop = FALSE]
    
    f_list[[kk]] = Z_tilde%*%(W_list[[kk]])
  }
  
  d_list[1]  = 0
  
  d_hat  =  d_list[Ktilde]
  W_hat  =  W_list[[Ktilde]]
  f_hat  =  f_list[[Ktilde]]
  
  return(list(d_hat   = d_hat,
              W_hat   = W_hat,
              f_hat   = f_hat,
              d_list  = d_list,
              W_list  = W_list,
              f_list  = f_list))
}


est.d.Wf.nPQ = function(Z, Ktilde = 10){
  n = dim(Z)[1];d1 = dim(Z)[2];d2 = dim(Z)[3];
  
  Z_tilde = Vec.tensor(Z)
  
  M = 0
  W_list = f_list = list()
  d_list = vector()
  dmax = d1*d2
  dstar = max(d1,d2)
  for (kk in 1:Ktilde){
    
    S_ztilde_k = sigmak(t(Z_tilde),as.matrix(colMeans(Z_tilde)),n = n, k= kk)
    
    M = M + S_ztilde_k%*%t(S_ztilde_k)
    
    ev_M = eigen(M)
    
    evalues = ev_M$values
    
    d_list[kk] =  max(which.max(evalues[1:(dmax-1)]/evalues[2:(dmax)]),dstar)
    
    W_list[[kk]] = ev_M$vectors[,1:d_hat, drop = FALSE]
    
    f_list[[kk]] = Z_tilde%*%(W_list[[kk]])
  }
  
  d_list[1]  = 0
  
  
  d_hat  =  d_list[Ktilde]
  W_hat  =  W_list[[Ktilde]]
  f_hat  =  f_list[[Ktilde]]
  
  
  return(list(d_hat   = d_hat,
              W_hat   = W_hat,
              f_hat   = f_hat,
              d_list  = d_list,
              W_list  = W_list,
              f_list  = f_list))
}


est.Wf = function(Y,P,Q,d,Ktilde = 10, thresh = FALSE, delta = NULL){
  n = dim(Y)[1];p = dim(Y)[2];q = dim(Y)[3];
  d1 = NCOL(P);d2 = NCOL(Q);
  
  Z = array(NA,dim = c(n,d1,d2))
  for (tt in 1:n) {
    Z[tt,,] = t(P)%*%Y[tt,,]%*%Q
  }
  
  Z_tilde = matrix(NA,n,d1*d2)
  for (tt in 1:n) {
    Z_tilde[tt,] = as.vector(Z[tt,,])
  }
  
  M = 0
  W_list = f_list = list()
  
  for (kk in 1:Ktilde){
    if(thresh){
      Y_2d <- Vec.tensor(Y)
      S_ytilde_k <- sigmak(t(Y_2d), as.matrix(colMeans(Y_2d)), n = n, k = kk)
      S_ytilde_k <- thresh_C(S_ytilde_k, delta)
      S_ztilde_k <- MatMult(MatMult(t(Q) %x% t(P), S_ytilde_k), Q %x% P)
    }
    else{S_ztilde_k = sigmak(t(Z_tilde),as.matrix(colMeans(Z_tilde)),n = n, k= kk)}
    
    M = M + S_ztilde_k%*%t(S_ztilde_k)
    
    ev_M = eigen(M)
    
    W_list[[kk]] = ev_M$vectors[,1:d, drop = FALSE]
    f_list[[kk]] = Z_tilde%*%(W_list[[kk]])
  }
  
  W_hat  =  W_list[[Ktilde]]
  f_hat  =  f_list[[Ktilde]]
  return(list(W_hat   = W_hat,
              f_hat   = f_hat,
              W_list  = W_list,
              f_list  = f_list))
}


est.UV.JAD = function(W,d1,d2,d){
  
  W_tilde_tol = array(NA,dim= c(d,d1,d2))
  
  for (jj in 1:d){
    W_tilde_i = matrix(NA,d1,d2)
    for (pp in 1:d2){
      for (mm in 1:d1) {
        W_tilde_i[mm,pp] <- W[mm + (pp - 1)*d1,jj]
      }
    }
    W_tilde_tol[jj,,] = W_tilde_i
  }
  
  P_tol = vector()
  for (ss in 1:d) {
    for(rr in ss:d){
      if(ss == rr){
        P_rs = minor_P(W_tilde_tol[rr,,],W_tilde_tol[ss,,],d1,d2)
      }else{
        P_rs = minor_P(W_tilde_tol[rr,,],W_tilde_tol[ss,,],d1,d2)
      }
      P_tol = cbind(P_tol,P_rs)
    }
  }
  
  M_tol    = svd(P_tol)$v
  dt       = NCOL(M_tol)
  M        = M_tol[,(dt-d+1):dt]
  M_tensor = array(NA,dim = c(d,d,d))
  eigen_gap = vector()
  for (ii in 1:d) {
    M_tensor[,,ii] = Vech2Mat_new(M[,ii], d)
    eigen_gap[ii] = min(abs(eigen(M_tensor[,,ii])$values))
  }
  
  
  ##construct H_star
  Ms = M_tensor[,,which.max(eigen_gap)]
  
  HH0 = vector()
  HH1 = vector()
  HH2 = vector()
  
  for (ii in 1:d) {
    HH0 = cbind(HH0,c(M_tensor[,,ii]))
    HH1 = cbind(HH1,c(solve(Ms)%*%M_tensor[,,ii]))
    HH2 = cbind(HH2,c(M_tensor[,,ii]%*%solve(Ms)))
  }
  
  PP = (t(HH0)%*%HH2)%*%solve(t(HH1)%*%HH2 + t(HH2)%*%HH1)%*%(t(HH2)%*%HH0)
  
  EVD = eigen(PP)
  if( min(EVD$values) < 0){
    R_NOJD = EVD$vectors
  }else{
    R_NOJD = sqrt(2)/2*(EVD$vectors)%*%diag((EVD$values)^{-1/2})%*%t(EVD$vectors)
  }
  
  M_1 = M%*%R_NOJD
  
  M_tensor_1 = array(NA,dim = c(d,d,d))
  for (ii in 1:d) {
    M_tensor_1[,,ii] = Vech2Mat_new(M_1[,ii],d)
    
  }
  
  H = jointDiag::ffdiag(M_tensor_1)$B ## jointDiag::ffdiag
  
  Theta = apply(MASS::ginv(H), 2, l2s)
  
  Wt = W%*%Theta
  
  U = matrix(NA,d1,d)
  V = matrix(NA,d2,d)
  
  for (jj in 1:d){
    Wt_tilde_i = matrix(NA,d1,d2)
    for (pp in 1:d2){
      for (mm in 1:d1) {
        Wt_tilde_i[mm,pp] <- Wt[mm + (pp - 1)*d1,jj]
      }
    }
    svdi = svd(Wt_tilde_i)
    U[,jj] = svdi$u[,1]
    V[,jj] = svdi$v[,1]
  }
  
  return(list(U = U,V = V,Theta = Theta))
  
}

est.UV.EVD = function(W,d1,d2,d){
  
  W_tilde_tol = array(NA,dim= c(d,d1,d2))
  
  for (jj in 1:d){
    W_tilde_i = matrix(NA,d1,d2)
    for (pp in 1:d2){
      for (mm in 1:d1) {
        W_tilde_i[mm,pp] <- W[mm + (pp - 1)*d1,jj]
      }
    }
    W_tilde_tol[jj,,] = W_tilde_i
  }
  
  P_tol = vector()
  for (ss in 1:d) {
    for(rr in ss:d){
      if(ss == rr){
        P_rs = minor_P(W_tilde_tol[rr,,],W_tilde_tol[ss,,],d1,d2)
      }else{
        P_rs = minor_P(W_tilde_tol[rr,,],W_tilde_tol[ss,,],d1,d2)
      }
      P_tol = cbind(P_tol,P_rs)
    }
  }
  Theta = 1 + 1i
  M_tol    = svd(P_tol)$v
  dt       = NCOL(M_tol)
  M_org        = M_tol[,(dt-d+1):dt]
  
  
  M  = M_org
  
  M_tensor = array(NA,dim = c(d,d,d))
  
  for (ii in 1:d) {
    M_tensor[,,ii] = Vech2Mat_new(M[,ii],d)
  }
  
  L0 = L1 = 0
  for (tt in 1:(d-1)) {
    L0 = L0 + M_tensor[,,tt]
    L1 = L1 + M_tensor[,,tt + 1]
  }
  
  L0 = L0/(d-1)
  L1 = L1/(d-1)
  
  theta.l = eigen(solve(L1)%*%L0)$vectors
  
  Theta = apply(L0%*%theta.l,2,l2s)
  
  if(is.complex(Theta)){
    Theta = Complex2Real(Theta)
  }
  
  Wt = W%*%Theta
  
  U = matrix(NA,d1,d)
  V = matrix(NA,d2,d)
  
  for (jj in 1:d){
    Wt_tilde_i = matrix(NA,d1,d2)
    for (pp in 1:d2){
      for (mm in 1:d1) {
        Wt_tilde_i[mm,pp] <- Wt[mm + (pp - 1)*d1,jj]
      }
    }
    svdi = svd(Wt_tilde_i)
    U[,jj] = svdi$u[,1]
    V[,jj] = svdi$v[,1]
  }
  
  return(list(U = U,V = V,Theta = Theta))
  
}

Sigma_Ycheck <- function(Y, k){
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  q <- dim(Y)[3]
  sigmaYk <- matrix(0, nrow = p^2 * q, ncol = q)
  Y_mean <- apply(Y, c(2, 3), mean)
  for (t in (k + 1):n) {
    A <- Y[t, , ] - Y_mean # Y_t - Y_bar (p x q)
    B <- Y[t - k, , ] - Y_mean # Y_{t-k} - Y_bar (p x q)
    
    # vec(B): 将 B 转换为列向量
    B_vec <- as.vector(B)
    
    # Kronecker product
    kron_prod <- A %x% B_vec #  (p*q) x (p*q)
    sigmaYk <- sigmaYk + kron_prod
  }
  return(sigmaYk/(n-k))
}

adjust_sign <- function(column) {
  first_nonzero_idx <- which(column != 0)[1]
  if (!is.na(first_nonzero_idx)) {
    return(sign(column[first_nonzero_idx]))
    } else {
      return(1)
    }
}