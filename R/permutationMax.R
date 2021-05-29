#' @importFrom stats acf ar 
#' @importFrom tseries garch
#' @useDynLib HDTSA
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @import Rcpp 
#' @import RcppEigen


## Output1: $NoGroups -- No. of groups with at least two components series
## Output2: $Nos_of_Members -- number of members in each of groups
## 				     with at least two members
## Output3: $Groups -- indices of components in each of groups with at least two members
## Output4: $maxcorr -- maximum correlation (over lags) of \eqn{p(p-1)/2} pairs in descending order
## Output5: $corrRatio -- ratios of succesive values from $maxCorr
## Output6: $NoConnectedPairs -- No. connected Pairs
#
#

permutationMax <- function(X, isvol=FALSE, m=NULL) {
  #
  # X: nxp data matrix
  # m: maximum lag used in calculating cross correlation coefficients

  ## Step 0: prewhiten each columns of X
  if(is.null(m)) m <- 10*log10(n/p)
  
  p=ncol(X) 
  n=nrow(X)
  if(!isvol)
  {
    R <- 5
    arOrder <- rep(0, p)
    for(j in 1:p) { 
      t=ar(X[,j], order.max = R)
      X[,j] <- t$resid; arOrder[j] <- t$order
      }
    j <- max(arOrder)
    X <- X[(j+1):n, ]
  }
  if(isvol)
  {
    ## Step 0: prewhiten each columns of X
    nanum=rep(0,p)
    for(j in 1:p) {options( warn = -1 )
                   t=tseries::garch(X[,j], order = c(1,1),trace=FALSE)
                    X[,j]=t$residuals
                    a=X[,j]
                    nanum[j]=length(a[is.na(X[,j])])

    }

    X=X[(max(nanum)+1):n,]
  }
  ## Step 1: calculate max_k |\rho(k)| for each pair components

  rho=acf(X,lag.max=m, plot=F) # rho$acf is an (m+1)xpxp array

  p0=p*(p-1)/2 # total number of pairs of component series
  M=vector(mode="numeric",length=p0) # max correlations between i-th and j-th component
                           # over lags between -m to m, for 1 <= j < i <= p
  for(i in 2:p) { for(j in 1:(i-1))
  	M[(i-2)*(i-1)/2+j]=max(abs(rho$acf[,i,j]), abs(rho$acf[,j,i]))
  	}
      # For a pxp matrix,  stack rows below the main diagoal together,
      # the (i,j)-th element, for i>j, is in the position (i-2)*(i-1)/2+j
  # cat("STEP1","\n")

  ## Step 2: sorting p0 maximum correlation in descending order,
  ## find the ratio estimator for r
  Ms=sort.int(M, decreasing=T, index.return=T)
      # Ms$x are sorted correlation, Ms$ix are the corresponding indices in M
  p1=as.integer(p0*0.75)
  ratio=Ms$x[1:p1]/Ms$x[2:(p1+1)]
  max0=max(ratio)
  n=1:p1
  r=n[ratio[n]==max0]; j=length(r); r=r[j]
  # cat("STEP2","\n")

  ## Step 3: find the pairs corresponding to the r maximum max_k |\rho(k)|
  h=mat.or.vec(p,1)
  for(i in 2:p) h[i]=(i-2)*(i-1)/2
  Inx=mat.or.vec(p,p); I=2:p
  for(k in 1:r) { q=I[(Ms$ix[k]-h[I])>0]; s=length(q); i=q[s]; j=Ms$ix[k]-h[i]; Inx[i,j]=1}
  # Now the entrices of Inx equal 1 are the positions with (i,j) connected,
  # and all other entrices are 0
  # cat("STEP3","\n")

  ## Step 4: picking up the grouping from each columns of Inx, mark column with Index=1
  ##         with a group with at least two members, and Index=0 otherwise
  G=mat.or.vec(p,p-1);
  Index=rep(0,p-1)
  N=rep(0,p-1)
  # G[,j] records the components (from j-th column of Inx) to be grouped together with j
  for(j in 1:(p-1)) { k=1
  	for(i in (j+1):p) if(Inx[i,j]>0) { k=k+1; G[k,j]=i}
  	if(k>1) { G[1,j]=j; Index[j]=1; N[j]=k }
  }
  # cat("STEP4","\n")

  ## Step 5: combining together any two groups obtained in Step 4 sharing
  ##         the same component
  check=1
  while(check>0) {  check=0
  for(i in 1:(p-2)) { if(Index[i]==0) next
          for(j in (i+1):(p-1)) { if(Index[j]==0) next
                  a=G[,i][G[,i]>0]; b=G[,j][G[,j]>0]
                  c=c(a, b); d=length(c[duplicated(c)])
                  if(d>0) {  # there are duplicated elements in a & b
                             a=unique(c) # picking up different elements from G[,i] & G[,j]
                             G[,i]=0; G[1:length(a),i]=sort(a); Index[j]=0
  			   N[i]=length(a); N[j]=0;
  			   check=1;
  		}
          }
  }
  # cat("d=", d, "\n")
  }
  # cat("STEP5","\n")


  ## Step 6: Output
  K=length(Index[Index==1])
  if(K==0) stop("All component series are linearly independent")
  Group=mat.or.vec(p,K+1)
  k=1
  for(j in 1:(p-1)) { if(Index[j]==1) { Group[,k]=G[,j]; k=k+1} }
  cat("\n"); cat("No of groups with more than one members:", K, "\n")
  cat("Nos of members in those groups:", N[N>0], "\n")
  cat("No of connected pairs:", r, "\n")
  #cat("Prewhited data are saved in the file Xpre.dat","\n\n")
  #write.table(X, "Xpre.dat", row.names=F, col.names=F)
  output=list(NoGroups=K, Nos_of_Members=N[N>0], Groups=Group[,1:K], maxcorr=Ms$x, corrRatio=ratio,NoConnectedPairs=r,Xpre=X)
  return(output)
}

