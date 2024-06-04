#' @title Estimation of matrix CP-factor model
#' @description \code{CP_MTS()} deals with CP-decomposition for high-dimensional
#'  matrix time series proposed in Chang et al. (2023):\deqn{{\bf{Y}}_t = {\bf A \bf X}_t{\bf B}^{'} +
#' {\boldsymbol{\epsilon}}_t, } where \eqn{{\bf X}_t = diag(x_{t,1},\ldots,x_{t,d})} is an \eqn{d \times d}
#' latent process, \eqn{{\bf A}} and \eqn{{\bf B}} are , respectively, \eqn{p
#' \times d} and \eqn{q \times d} unknown constant matrix, and \eqn{ {\boldsymbol{\epsilon}}_t }
#'  is a \eqn{p \times q} matrix white noise process. This function aims to estimate the rank
#'  \eqn{d} and the coefficient matrices \eqn{{\bf A}} and \eqn{{\bf B}}.
#'
#' @param Y A \eqn{n \times p \times q} data array, where \eqn{n} is the sample size and \eqn{(p,q)}
#' is the dimension of \eqn{{\bf Y}_t}.
#' @param xi A \eqn{n \times 1} vector. If \code{NULL} (the default), then a PCA-based \eqn{\xi_{t}}
#' is used [See Section 5.1 in Chang et al. (2023)] to calculate the sample auto-covariance matrix
#' \eqn{\widehat{\bf \Sigma}_{\bf Y, \xi}(k)}.
#' @param Rank A list of the rank \eqn{d},\eqn{d_1} and \eqn{d_2}. Default to \code{NULL}.
#' @param lag.k Integer. Time lag \eqn{K} is only used in \code{CP.Refined} and \code{CP.Unified} to
#' calculate the nonnegative definte matrices \eqn{\widehat{\mathbf{M}}_1} and
#' \eqn{\widehat{\mathbf{M}}_2}: \deqn{\widehat{\mathbf{M}}_1\ =\
#'   \sum_{k=1}^{K}\widehat{\mathbf{\Sigma}}_{\bf Y, \xi}(k)\widehat{\mathbf{\Sigma}}_{\bf Y, \xi}(k)',
#'   }, \deqn{\widehat{\mathbf{M}}_2\ =\
#'   \sum_{k=1}^{K}\widehat{\mathbf{\Sigma}}_{\bf Y, \xi}(k)'\widehat{\mathbf{\Sigma}}_{\bf Y, \xi}(k),
#'   }
#'   where \eqn{\widehat{\mathbf{\Sigma}}_{\bf Y, \xi}(k)} is the sample auto-covariance of
#'   \eqn{ {\bf Y}_t} and \eqn{\xi_t} at lag \eqn{k}.
#' @param lag.ktilde Integer. Time lag \eqn{\tilde K} is only used in \code{CP.Unified} to calulate the
#' nonnegative definte matrix \eqn{\widehat{\mathbf{M}}}: \deqn{\widehat{\mathbf{M}} \ =\
#'   \sum_{k=1}^{\tilde K}\widehat{\mathbf{\Sigma}}_{\tilde{\bf Z}}(k)\widehat{\mathbf{\Sigma}}_{\tilde{\bf Z}}(k)'.
#'   }
#' @param method Method to use: \code{CP.Direct} and \code{CP.Refined}, Chang et al.(2023)'s direct and refined estimators;
#'  \code{CP.Unified}, Chang et al.(2024+)'s unified estimation procedure.
#'
#' @return An object of class "mtscp" is a list containing the following
#'   components:
#'   \item{A}{The estimated \eqn{p \times d} left loading matrix \eqn{\widehat{\bf A}}.}
#'   \item{B}{The estimated \eqn{q \times d} right loading matrix \eqn{\widehat{\bf B}}.}
#'   \item{f}{The estimated latent process \eqn{(\hat{x}_{1,t},\ldots,\hat{x}_{d,t})}.}
#'   \item{Rank}{The estimated rank \eqn{(\hat{d}_1,\hat{d}_2,\hat{d})} of the matrix CP-factor model.}
#'
#'
#' @references
#'   Chang, J., He, J., Yang, L. and Yao, Q.(2023). \emph{Modelling matrix time series via a tensor CP-decomposition}.
#'   Journal of the Royal Statistical Society Series B: Statistical Methodology, Vol. 85(1), pp.127--148.
#'   
#'   Chang, J., Du, Y., Huang, G. and Yao, Q.(2024+). \emph{On the Identification and Unified Estimation
#'   Procedure for the Matrix CP-factor Model}, Working paper.
#'
#' @examples
#' p = 10
#' q = 10
#' n = 400
#' d = d1 = d2 = 3
#' data <- DGP.CP(n,p,q,d,d1,d2)
#' Y = data$Y
#' res1 <- CP_MTS(Y,method = "CP.Direct")
#' res2 <- CP_MTS(Y,method = "CP.Refined")
#' res3 <- CP_MTS(Y,method = "CP.Unified")
#' @export
#' @useDynLib HDTSA
#' @importFrom stats arima.sim rnorm runif

CP_MTS = function(Y,xi = NULL, Rank = NULL, lag.k = 15, lag.ktilde =  10, method = c("CP.Direct","CP.Refined","CP.Unified")){
  n = dim(Y)[1]; p = dim(Y)[2]; q = dim(Y)[3];
  if(is.null(xi)){
    xi = est.xi(Y)
  }
  if(method == "CP.Direct"){
    S_yxi_1 = Autocov_xi_Y(Y,xi,lag.k = 1)
    S_yxi_2 = Autocov_xi_Y(Y,xi,lag.k = 2)
    if(p > q){
      ##(1) estimation of d
      K1 = t(S_yxi_1)%*%S_yxi_1
      eg1 = eigen(K1)
      w = eg1$values
      ww = w[-1]/w[-length(w)]
      d = which(ww==min(ww[1:floor(0.75*q)]))

      if (d > 1){
        K1 = eg1$vectors[,1:d]%*%diag(eg1$values[1:d])%*%t(eg1$vectors[,1:d]);
      }else{
        K1 = eg1$vectors[,1]%*%diag(eg1$values[1],1)%*%t(eg1$vectors[,1]);
      }
        K2 = t(S_yxi_1)%*%S_yxi_2;

      ##(2) estimation of A and B
      Geg = geigen::geigen(K2,K1);
      evalues = Geg$values[which(Mod(Geg$values)<=10^5&Geg$values!=0)]

      Bl = Geg$vectors[,which(Geg$values %in% evalues)]
      A = apply(S_yxi_1%*%Bl,2,l2s)
      Al = t(MASS::ginv(A))
      B = apply(t(S_yxi_1)%*%Al,2,l2s)
    }else{

      ##(1) estimation of d
      K1 = S_yxi_1%*%t(S_yxi_1)
      eg1 = eigen(K1)
      w = eg1$values
      ww = w[-1]/w[-length(w)]
      d = which(ww==min(ww[1:floor(0.75*p)]))

      if (d > 1){
        K1 = eg1$vectors[,1:d]%*%diag(eg1$values[1:d])%*%t(eg1$vectors[,1:d]);
      }else{
        K1 = eg1$vectors[,1]%*%diag(eg1$values[1],1)%*%t(eg1$vectors[,1]);
      }
      K2 = S_yxi_1%*%t(S_yxi_2);
      ##(2) estimation of A and B
      Geg = geigen::geigen(K2,K1);
      evalues = Geg$values[which(Mod(Geg$values)<=10^5&Geg$values!=0)]
      Al = Geg$vectors[,which(Geg$values %in% evalues)]
      B = apply(t(S_yxi_1)%*%Al,2,l2s)
      Bl = t(MASS::ginv(B))
      A = apply(S_yxi_1%*%Bl,2,l2s)
    }
    ##(3) estimation of Xt
    H = matrix(NA,p*q,d)
    for (ii in 1:d) {
      H[,ii] = B[,ii] %x% A[,ii]
    }
    f = Vec.tensor(Y)%*%H%*%MASS::ginv(t(H)%*%H)

    if(is.complex(A) == T || is.complex(B) == T ){
      A = Complex2Real(A)
      B = Complex2Real(B)
      f = Complex2Real(f)
    }
    METHOD <- c("Estimation of matrix CP-factor model",paste("Method:",method))
    names(d) <- "The estimated rank of the matrix CP-factor model"
    con = structure(list(A = A,B = B,f = f,Rank = d, method = METHOD),
                    class = "mtscp")
    return(con)
  }
  if(method == "CP.Refined"){
    ##(1) estimation of P,Q and d
    M1 = M2 = 0
    dmax = round(min(p,q)*0.75)

    for (kk in 1:lag.k){
      S_yxi_k = Autocov_xi_Y(Y,xi,lag.k = kk)
      M1 = M1 + S_yxi_k%*%t(S_yxi_k)
      M2 = M2 + t(S_yxi_k)%*%S_yxi_k
    }
    ev_M1 = eigen(M1)
    ev_M2 = eigen(M2)

    d1 =  which.max(ev_M1$values[1:dmax]/ev_M1$values[2:(dmax+1)])
    d2 =  which.max(ev_M2$values[1:dmax]/ev_M2$values[2:(dmax+1)])

    d   = ifelse(p>q,d1,d2)

    if(!is.null(Rank)){
      d = Rank
    }
      P = ev_M1$vectors[,1:d]
      Q = ev_M2$vectors[,1:d]

    if(d == 1){
      A = as.matrix(P)
      B = as.matrix(Q)

      f = vector()
      for (tt in 1:n) {
        f[tt] = t(A)%*%Y[tt,,]%*%B
      }
      f = as.matrix(f)

    }else{
      ##(2) estimation of U and V
      Z = array(NA,dim = c(n,d,d))
      for (tt in 1:n) {
        Z[tt,,] = t(P)%*%Y[tt,,]%*%Q
      }
      xi =  est.xi(Z)
      S_Zxi_1 = Autocov_xi_Y(Z,xi,lag.k = 1)
      S_Zxi_2 = Autocov_xi_Y(Z,xi,lag.k = 2)

      vl = eigen(MASS::ginv(t(S_Zxi_1)%*%S_Zxi_1)%*%t(S_Zxi_1)%*%S_Zxi_2)$vectors ##MASS
      ul = eigen(MASS::ginv(S_Zxi_1%*%t(S_Zxi_1))%*%S_Zxi_1%*%t(S_Zxi_2))$vectors

      U = apply(S_Zxi_1%*%vl,2,l2s)
      V = apply(t(S_Zxi_1)%*%ul,2,l2s)

      ##(3) estimation of A and B
      A = P%*%U
      B = Q%*%V

      ##(4) estimation of Xt
      W = matrix(NA,d^2,d)

      for (ii in 1:d) {
        W[,ii] = V[,ii]%x%U[,ii]
      }

      f = Vec.tensor(Z)%*%W%*%solve(t(W)%*%W)

      if(is.complex(A) == T || is.complex(B) == T ){
        A = Complex2Real(A)
        B = Complex2Real(B)
        f = Complex2Real(f)
      }
    }
      METHOD <- c("Estimation of matrix CP-factor model",paste("Method:",method))
      names(d) <- "The estimated rank of the matrix CP-factor model"
    con = structure(list(A = A,B = B,f = f,Rank = d, method = METHOD),
                    class = "mtscp")
    return(con)
  }
  if(method == "CP.Unified"){
    ##(1) estimation of P,Q and d1,d2
    if(is.null(Rank)){

      PQ_hat_tol = est.d1d2.PQ(Y,xi,K = lag.k)

      d1  = PQ_hat_tol$d1_hat
      d2  = PQ_hat_tol$d2_hat
      d   = NULL
      P   = PQ_hat_tol$P_hat
      Q   = PQ_hat_tol$Q_hat

      if(d1 == 1 || d2 == 1){d = d1*d2}

    }else{

      d   = Rank$d
      d1  = Rank$d1
      d2  = Rank$d2

      PQ_hat_tol = est.PQ(Y,xi,d1,d2,K = lag.k)

      P   = PQ_hat_tol$P_hat
      Q   = PQ_hat_tol$Q_hat

    }
    ##(2) estimation of W* = (v1*u1,v2*u2,...,vd*ud)H = WH
    if(d1 == 1 & d2 == 1){
      d = 1
      f = vector()
      for (tt in 1:n) {
        f[tt] = t(P)%*%Y[tt,,]%*%Q
      }
      A = P
      B = Q

      # con = list(A = A,B = B,f = f,Rank = list(d=d,d1=d1,d2=d2))
      METHOD <- c("Estimation of matrix CP-factor model",paste("Method:",method))
      rank <- list(d=d,d1=d1,d2=d2)
      names(rank) <- c("The estimated rank d of the matrix CP-factor model",
                       "The estimated rank d1 of the matrix CP-factor model",
                       "The estimated rankd d2 of the matrix CP-factor model")
      con = structure(list(A = A,B = B,f = f,Rank = rank, method = METHOD),
                      class = "mtscp")

      return(con)
    }else{

      if(is.null(d)){
        W_hat_tol  =  est.d.Wf(Y,P,Q,Ktilde =  lag.ktilde)
        d          =  W_hat_tol$d_hat
        W          =  W_hat_tol$W_hat
        f          =  W_hat_tol$f_hat
      }else{
        d          =  d
        W_hat_tol  =  est.Wf(Y,P,Q,d,Ktilde = lag.ktilde)
        W          =  W_hat_tol$W_hat
        f          =  W_hat_tol$f_hat
      }

      ##(3) estimation of U and V
      if(d1 == 1 || d2 == 1){
        Theta = NULL
        if(d1 == 1){
          U = 1;
          V = W;
          warning("d1 equal to 1, V cannot be identified uniquely!")
        }
        if(d2 == 1){
          U = W;
          V = 1;
          warning("d2 equal to 1, U cannot be identified uniquely!")
        }
        if(d1 == 1 & d2 == 1){
          U = 1;
          V = 1;
        }
        U = as.matrix(U)
        V = as.matrix(V)
      }else{

        UV_hat_tol = est.UV.JAD(W,d1,d2,d)

        U          = UV_hat_tol$U
        V          = UV_hat_tol$V
        Theta      = UV_hat_tol$Theta
      }


      ##(4) estimation of A and B
      A = P%*%U
      B = Q%*%V
      METHOD <- c("Estimation of matrix CP-factor model",paste("Method:",method))
      rank <- list(d=d,d1=d1,d2=d2)
      names(rank) <- c("The estimated rank d of the matrix CP-factor model",
                       "The estimated rank d1 of the matrix CP-factor model",
                       "The estimated rankd d2 of the matrix CP-factor model")
      con = structure(list(A = A,B = B,f = f,Rank = rank, method = METHOD),
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


#' @title Data generate process of matrix CP-factor model
#' @description \code{DGP.CP()} function generate the matrix time series described in Chang et al. (2023):\deqn{{\bf{Y}}_t = {\bf A \bf X}_t{\bf B}^{'} +
#' {\boldsymbol{\epsilon}}_t, } where \eqn{{\bf X}_t = diag(x_{t,1},\ldots,x_{t,d})} is an \eqn{d \times d}
#' latent process, \eqn{{\bf A}} and \eqn{{\bf B}} are , respectively, \eqn{p
#' \times d} and \eqn{q \times d} unknown constant matrix, and \eqn{ {\boldsymbol{\epsilon}}_t }
#'  is a \eqn{p \times q} matrix white noise process.
#'
#' @param n  Integer. Sample size of \eqn{\bf Y_t}, \eqn{t=1,\ldots,n}.
#' @param p  Integer. Number of rows of \eqn{\bf Y_t}.
#' @param q  Integer. Number of columns of \eqn{\bf Y_t}.
#' @param d1  Integer. Rank of \eqn{\bf A}.
#' @param d2  Integer. Rank of \eqn{\bf B}.
#' @param d  Integer. Number of columns of \eqn{\bf A} and \eqn{\bf B}.
#'
#' @seealso \code{\link{CP_MTS}}.
#' @return A list containing the following
#'   components:
#'   \item{Y}{A \eqn{n \times p \times q} data array of \eqn{\bf Y_t}.}
#'   \item{S}{A \eqn{n \times p \times q} data array of \eqn{\bf S_t = \bf A \bf X_t \bf B'}.}
#'   \item{A}{A \eqn{p \times d} coefficient matrix.}
#'   \item{B}{A \eqn{q \times d} coefficient matrix.}
#'   \item{X}{A \eqn{n \times d \times d} data array of \eqn{\bf X_t}.}
#'   \item{P}{A \eqn{p \times d_1} orthogonal matrix such that \eqn{\bf A = \bf P \bf U}.}
#'   \item{Q}{A \eqn{q \times d_2} orthogonal matrix such that \eqn{\bf B = \bf Q \bf V}.}
#'   \item{U}{A \eqn{d_1 \times d} matrix such that \eqn{\bf A = \bf P \bf U}.}
#'   \item{V}{A \eqn{d_2 \times d} matrix such that \eqn{\bf B = \bf Q \bf V}.}
#'   \item{W}{A \eqn{d_1 d_2 \times d} matrix such that \eqn{\bf W = (\bf v_1 \otimes \bf u_1,\ldots,\bf v_d \otimes \bf u_d)}.}
#'   \item{Ws}{A \eqn{d_1 d_2 \times d} matrix. An orthogonal basis of \eqn{\bf W}.}
#'   \item{Xmat}{A \eqn{n \times d} data matrix of \eqn{diag(\bf X_t)}.}
#'   \item{Smat}{A \eqn{n \times pq} data matrix of \eqn{vec(\bf S_t)}.}
#' @references
#'   Chang, J., He, J., Yang, L. and Yao, Q.(2023). \emph{Modelling matrix time series via a tensor CP-decomposition}.
#'   Journal of the Royal Statistical Society Series B: Statistical Methodology, Vol. 85(1), pp.127--148.
#'
#' @examples
#' p = 10
#' q = 10
#' n = 400
#' d = d1 = d2 = 3
#' data <- DGP.CP(n,p,q,d,d1,d2)
#' Y = data$Y
#' @export
DGP.CP = function(n,p,q,d1,d2,d){

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
              S      =  S,
              A      =  A_s,
              B      =  B_s,
              X      =  X,
              P      =  P,
              Q      =  Q,
              U      =  U,
              V      =  V,
              W      =  W,
              Ws     =  W_star,
              Xmat    =  X_m,
              Smat    =  S_m,
              Input  =  Input
  ))

}


Autocov_xi_Y = function(Y,xi,lag.k = k){

  n = dim(Y)[1]
  k = lag.k

  Y_mean = 0
  xi_mean = 0
  for (ii in 1:n) {
    Y_mean = Y_mean + Y[ii,,]
    xi_mean = xi_mean + xi[ii]
  }
  Y_mean   = Y_mean/n
  xi_mean = xi_mean/n

  Sigma_Y_xi_k = 0
  for (ii in (k+1):n) {
    Sigma_Y_xi_k = Sigma_Y_xi_k + (Y[ii,,] - Y_mean)*(xi[ii-k] - xi_mean)
  }
  return(Sigma_Y_xi_k/(n-k))
}

est.xi  = function(Y, thresh_per = 0.99){

  n = dim(Y)[1];p = dim(Y)[2];q = dim(Y)[3];

  xi.mat = Vec.tensor(Y)

  xi.mat = scale(xi.mat,scale = F)

  if(n > p*q){
    eig_xi.mat = eigen(MatMult(t(xi.mat),xi.mat))
    cfr =  cumsum(eig_xi.mat$values)/sum(eig_xi.mat$values)
    d_hat = min(which(cfr > thresh_per))
    w_hat = eig_xi.mat$vectors[,1:d_hat]

    xi.f = xi.mat%*%w_hat

    xi   = rowMeans(xi.f)

  }else{
    eig_xi.mat = eigen(MatMult(xi.mat,t(xi.mat)))
    cfr = cumsum(eig_xi.mat$values)/sum(eig_xi.mat$values)
    d_hat = min(which(cfr > thresh_per))
    xi.f1 = as.matrix(eig_xi.mat$vectors[,1:d_hat])

    weight = sqrt(eig_xi.mat$values[1:d_hat])

    xi.f1 = xi.f1%*%diag(weight)

    xi = rowMeans(xi.f1)

  }
  return(xi)
}

est.d1d2.PQ = function(Y,xi,K = 10){
  n = dim(Y)[1];p = dim(Y)[2];q = dim(Y)[3];

  d2_list = d1_list =vector()
  M1 = M2 = 0
  dmax = round(min(p,q)*0.75)
  P_list = Q_list = list()
  for (kk in 1:K){

    S_yxi_k = Autocov_xi_Y(Y,xi,lag.k = kk)

    M1 = M1 + S_yxi_k%*%t(S_yxi_k)
    M2 = M2 + t(S_yxi_k)%*%S_yxi_k

    ev_M1 = eigen(M1)
    ev_M2 = eigen(M2)

    d1_list[kk] =  which.max(ev_M1$values[1:dmax]/ev_M1$values[2:(dmax+1)])
    d2_list[kk] =  which.max(ev_M2$values[1:dmax]/ev_M2$values[2:(dmax+1)])

    P_list[[kk]] = ev_M1$vectors[,1:(d1_list[kk])]
    Q_list[[kk]] = ev_M2$vectors[,1:(d2_list[kk])]

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

est.PQ = function(Y,xi,d1,d2,K = 20){
  n = dim(Y)[1];p = dim(Y)[2];q = dim(Y)[3];

  M1 = M2 = 0
  dmax = round(min(p,q)*0.75)
  P_list = Q_list = list()
  for (kk in 1:K){

    S_yxi_k = Autocov_xi_Y(Y,xi,lag.k = kk)

    M1 = M1 + S_yxi_k%*%t(S_yxi_k)
    M2 = M2 + t(S_yxi_k)%*%S_yxi_k

    ev_M1 = eigen(M1)
    ev_M2 = eigen(M2)

    P_list[[kk]] = ev_M1$vectors[,1:(d1)]
    Q_list[[kk]] = ev_M2$vectors[,1:(d2)]

  }

  P_hat  =  P_list[[K]]
  Q_hat  =  Q_list[[K]]

  return(list(P_hat   = P_hat,
              Q_hat   = Q_hat,
              P_list  = P_list,
              Q_list  = Q_list))
}


est.d.Wf = function(Y,P,Q, Ktilde = 10){
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

    S_ztilde_k = sigmak(t(Z_tilde),as.matrix(colMeans(Z_tilde)),n = n, k= kk)

    M = M + S_ztilde_k%*%t(S_ztilde_k)

    ev_M = eigen(M)

    evalues = ev_M$values


    d_list[kk] =  max(which.max(evalues[1:(dmax-1)]/evalues[2:(dmax)]),dstar)

    W_list[[kk]] = ev_M$vectors[,1:(d_list[kk])]

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

    W_list[[kk]] = ev_M$vectors[,1:d_hat]

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


est.Wf = function(Y,P,Q,d,Ktilde = 10){
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
    S_ztilde_k = sigmak(t(Z_tilde),as.matrix(colMeans(Z_tilde)),n = n, k= kk)

    M = M + S_ztilde_k%*%t(S_ztilde_k)

    ev_M = eigen(M)

    W_list[[kk]] = ev_M$vectors[,1:d]
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
    M_tensor[,,ii] = Vech2Mat_new(M[,ii],d)
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
