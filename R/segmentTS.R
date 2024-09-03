#' @importFrom stats var
#' @importFrom clime clime cv.clime tracel2
#' @useDynLib HDTSA
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp


segmentTS <- function(Y, lag.k, 
                      thresh = FALSE, 
                      tuning.vec = NULL,
                      opt = 1,
                      control = list(),
                      K = 5)
{
  n=nrow(Y)
  p=ncol(Y)
  storage.mode(p)<-"integer"
  storage.mode(n)<-"integer"
  
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
      if(ev[i] <  1e-4)D[i, i] <- sqrt(log(p) / n)
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
  storage.mode(Y)<-"double"
  storage.mode(mean_y)<-"double"
  
  Wy=diag(rep(1,p))

  if(!thresh)
  {
    for(k in 1:lag.k) {
      Sigma_y<-sigmak(Y,mean_y,k,n)
      Wy=Wy+MatMult(Sigma_y,t(Sigma_y))
      #S=cov(t(Y[,1:(n-k)]),t(Y[,(1+k):n])); Wy=Wy+S%*%t(S)
    }
  }
  
  # Segment with thresholding
  
  if(thresh)
  {
    # Using the cross validation for thresholding
    
    for(k in 1:lag.k)
    {
      error=NULL
      storage.mode(k)<-"integer"
      if(is.null(tuning.vec))
      {
        tuning.vec=2
        deltafinal=tuning.vec
      }
      else{
        
        # To select proper threshold parameter
        if(length(tuning.vec)>1){
          for(v in 1:K)
          {
            sample1=sample(1:n,size=n/2)
            sample2=c(1:n)[-sample1]
            sampleY1=Y[,sample1]
            sampleY2=Y[,sample2]
            
            mean_y1<-as.matrix(rowMeans(sampleY1))
            mean_y2<-as.matrix(rowMeans(sampleY2))
            
            storage.mode(mean_y1)<-"double"
            storage.mode(sampleY1)<-"double"
            storage.mode(mean_y2)<-"double"
            storage.mode(sampleY2)<-"double"
            n1=ceiling(n/2)
            storage.mode(n1)<-"integer"
            storage.mode(p)<-"integer"
            
            errors=NULL
            
            for(d in 1:length(tuning.vec))
            {
              delta1=tuning.vec[d]
              storage.mode(delta1)<-"double"
              
              Sigma_y1 <- sigmak(sampleY1,mean_y1,k,n1)
              Sigma_y2 <- sigmak(sampleY2,mean_y2,k,n1)
              
              storage.mode(Sigma_y1)<-"double"
              storage.mode(Sigma_y2)<-"double"
              
              Sigma_ythres1 <- thresh_C(Sigma_y1, sampleY1, mean_y1, k, n1, p, delta1)
              
              errors=c(errors,(norm((Sigma_ythres1-Sigma_y2),type="F"))^2)
            }
            error=rbind(error,errors)
          }
          errormean=colMeans(error)
          d=which.min(errormean)
          deltafinal=tuning.vec[d]
        }
      }
      # Find the best tuning parameter
      
      #res<-.Fortran("segment",Y,mean_y,k,n,p,res=numeric(p^2))$res
      Sigma_y <- sigmak(Y,mean_y,k,n)
      storage.mode(Sigma_y)<-"double"
      storage.mode(deltafinal)<-"double"
      
      # Carry out the final thresholding
      
      Sigma_ynew <- thresh_C(Sigma_y, Y, mean_y, k, n, p, deltafinal)
      
      Wy=Wy+MatMult(Sigma_ynew,t(Sigma_ynew))
    }
  }
  # Segment without thresholding
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
