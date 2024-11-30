#' @importFrom stats var
#' @importFrom clime clime cv.clime tracel2
#' @useDynLib HDTSA
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp


segmentTS <- function(Y, lag.k, 
                      thresh = FALSE, 
                      delta = 2 * sqrt(log(ncol(Y)) / nrow(Y)),
                      opt = 1,
                      control = list()
                      )
{
  n=nrow(Y)
  p=ncol(Y)
  
  # Part I -- standardize Y  such that var(y_t)=I_p
  if(opt == 1)
  {
    M <- var(Y)
    t <- eigen(M, symmetric = T)
    ev <- t$values
    G <- as.matrix(t$vectors)
    D <- G * 0
    for (i in 1:p) {
      if (ev[i] > 1e-4)
        D[i, i] <- sqrt(1 / ev[i])
      else {
        D[i, i] <- 1 / sqrt(log(p) / n)
      }
    }
  }
  else if(opt == 2){
    #print('now use clime to calculate')
    con <- list(nlambda = 100, lambda.max = 0.8, lambda.min = ifelse(n>p,1e-4,1e-2),
                standardize = FALSE, linsolver = "primaldual")
    con[(namc <- names(control))] <- control
    
    
    M <- clime(Y, nlambda=con$nlambda, standardize = con$standardize,
               lambda.max = con$lambda.max, lambda.min = con$lambda.min,
               linsolver = con$linsolver)
    M <- clime::cv.clime(M, loss = "tracel2")
    #print(M$lambda)
    #print(M$lambdaopt)
    M <- clime(Y, standardize = con$standardize, lambda = M$lambdaopt)
    e <- unlist(M$Omega)
    names(e) <- NULL
    M <- matrix(e,nrow = p)
    t <- eigen(M,symmetric = T)
    G <- as.matrix(t$vectors)
    ev <- t$values
    D <- G * 0
    # square root of eigenvalues
    for(i in 1:p)
    {
      if(ev[i] <  1e-4 | ev[i] > 1e+4)D[i, i] <- sqrt(log(p) / n)
      else D[i, i]=sqrt(ev[i])
    }
  }

  # M1=var(y_t)^{-1/2}
  M1 <- MatMult(MatMult(G, D), t(G))
  Y1 <- MatMult(M1, t(Y))
  # Y is standardized now: var(y_t)=I_p
  Y <- Y1

  # Part II -- Apply the transformation to recover x_t

  
  mean_y<-as.matrix(rowMeans(Y))
  
  Wy=diag(rep(1,p))
  
  for (k in 1:lag.k) {
    Sigma_y<-sigmak(Y,mean_y,k,n)
    if(!thresh)
      Wy=Wy+MatMult(Sigma_y,t(Sigma_y))
    else {
      Sigma_ynew <- thresh_C(Sigma_y, delta)
      Wy=Wy+MatMult(Sigma_ynew,t(Sigma_ynew))
    }
  }
  
  t=eigen(Wy, symmetric=T)
  G=as.matrix(t$vectors)
  Y1=MatMult(t(G),Y)
  # segmented series
  Z=t(Y1)
  # transformation matrix x_t = B y_t, does not include permutation in Step 2
  B=MatMult(t(G),M1)
  Yt=list(B=B, Z=Z)
  return(Yt)
}
