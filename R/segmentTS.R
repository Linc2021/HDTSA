#' Segment Multivariate Time Series
#'
#' Calculate linear transformation of the \eqn{p}-variate time series \eqn{y_t} such that the transformed series \eqn{x_t=By_t} is segmented into several
#' lower-dimensional subseries, and those subseries are uncorrelated with
#' each other both contemporaneously and serially.
#'
#' When \eqn{p} is small, thresholding is not required. However, when \eqn{p} is large, it is necessary to use the thresholding method, see more information in Chang et al. (2014).
#'
#' @param Y   A data matrix with \eqn{n} rows and \eqn{p} columns, where \eqn{n} is the sample size and \eqn{p} is the dimension of the time series.
#' @param lag.k  A positive integer specified to calculate Wy. See (2.5) in Chang et al. (2018).
#' @param thresh   Logical. If \code{FALSE} (the default), no thresholding will be applied. If \code{TRUE}, a
#'                 thresholding method will be applied first to estimate Wy, see (3.4) and (3.5) in Chang et al. (2018).
#' @param tuning.vec  The value of thresholding parameter \eqn{\lambda}. The thresholding level is specified by
#'                  \deqn{ u = \lambda {(log p/n)^(1/2)}.}
#'                    Default value is 2. If \code{tuning.vec} is a vector, then a cross validation method proposed in Cai and Liu (2011) will be used
#'                    to choose the best tuning parameter.
#' @param isvol Logical. If \code{FALSE} (the default), it means that the process is stationary, otherwise it is Volatility process.
#' @param K   The number of folders used in the cross validation, the default is 5. It is required when \code{thresh} is \code{TRUE}.
#' @return An object of class "segmentTS" is a list containing the following components:
#' \item{B}{the \eqn{p} by \eqn{p} transformation matrix such that \eqn{x_t=By_t}}
#' \item{X}{the transformed series with \eqn{n} rows and \eqn{p} columns}
#' @note This is the first step to transform the time series. The second step is grouping the transformed time series, see \code{\link{permutationMax}}, \code{\link{permutationFDR}}.
#'
#'
#' @references Chang, J., Guo, B. & Yao, Q. (2018). \emph{Principal component analysis for second-order stationary vector time series}. The Annals of Statistics, Vol. 46, pp. 2094–2124.
#'
#'             Cai, T. and Liu, W. (2011). \emph{Adaptive thresholding for sparse covariance matrix estimation}.  Journal of the American Statistical Association 106: 672-684.
## @family aggregate functions
#' @seealso \code{\link{permutationMax}}, \code{\link{permutationFDR}}.
#' @importFrom stats var
#' @importFrom clime clime cv.clime
#' @useDynLib HDTSA
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @import Rcpp 
#' @import RcppEigen
#' @export






segmentTS <- function(Y, lag.k, isvol = FALSE, 
                      thresh = FALSE, 
                      tuning.vec = NULL, 
                      K = 5)
{
  n=nrow(Y)
  p=ncol(Y)
  storage.mode(p)<-"integer"
  storage.mode(n)<-"integer"
  
  # Part I -- standardize Y  such that var(y_t)=I_p
  if(p<n)
  {
    M <- var(Y)
    t <- eigen(M, symmetric = T)
    ev <- t$values
    G <- as.matrix(t$vectors)
    D <- G * 0
    for (i in 1:p) {
      if (ev[i] > 0)
        D[i, i] <- sqrt(1 / ev[i])
      else {
        D[i, i] <- 1 / sqrt(log(p) / n)
        #print("Data are degenerate")
        #quit("yes")
      }
    }
  }
  else{
    #print('now use clime to calculate')
    M <- clime(Y, nlambda=10, standardize = FALSE, linsolver = "simplex")
    M <- cv.clime(M, loss = "tracel2")
    #print(M$lambda)
    #print(M$lambdaopt)
    M <- clime(Y, lambda = M$lambdaopt)
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
      if(ev[i] < 0)D[i, i] <- sqrt(log(p) / n)
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
  
  if(!isvol){
    Wy=vol_wy(Y,drop(mean_y),lag.k,n,p)
  }
  else{
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
  }
  # Segment without thresholding
  t=eigen(Wy, symmetric=T)
  G=as.matrix(t$vectors)
  Y1=MatMult(t(G),Y)
  # segmented series
  X=t(Y1)
  # transformation matrix x_t = B y_t, does not include permutation in Step 2
  B=MatMult(t(G),M1)
  Yt=list(B=B, X=X)
  return(Yt)
}
