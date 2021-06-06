#' @importFrom stats acf ar pnorm var
#' @useDynLib HDTSA
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @import Rcpp 
#' 
permutationFDR <- function(X,prewhiten=TRUE, beta, m=NULL, verbose = FALSE) {
  # X: nxp data matrix
  # m: maximum lag used in multipe test for each pair components
  # beta: the error rate in FDR
  # R: upper bound of AR-order in prewhitenning

  ## Step 0: prewhiten each columns of X
  if(is.null(m)) m <- 10
  #if(is.null(beta)) beta <- 10^(-8)
  p=ncol(X)
  n=nrow(X)

  if(prewhiten)
  {
    R <- 5
    arOrder=rep(0, p)
    for(j in 1:p) { 
      t=ar(X[,j], order.max=R)
      X[,j]=t$resid
      arOrder[j]=t$order }
    j=max(arOrder)
    X=X[(j+1):n,]
    n=n-j
  }
  ## Step 1: for each 1\le j < i \le p, calculate P-value for the multiple test
  ## for H_0
  sqn=sqrt(n); m2=2*m+1
  rho=acf(X,lag.max=m, plot=F) # rho$acf is an (m+1)xpxp array
  p0=p*(p-1)/2 # total number of pairs of component series
  Pv=vector(mode="numeric",length=m2)
  M=vector(mode="numeric",length=p0) 
  # P-value for testing H_0 for i-th and j-th components
  for(i in 2:p) 
    for(j in 1:(i-1)) { 
      Pv[m2]=2*pnorm(-sqn*abs(rho$acf[1,i,j]))
      for(k in 1:(m-1)) { Pv[k]=2*pnorm(-sqn*abs(rho$acf[k+1,i,j]))
  			Pv[m+k]=2*pnorm(-sqn*abs(rho$acf[k+1,j,i]))
  	}
  	Pv=sort(Pv); k=1:m2; Pv=(Pv/k)*m2
  	M[(i-2)*(i-1)/2+j]=min(Pv) 
  	# P-value for multiple test for H_0 for (i,j) pairs
  }
      # For a pxp matrix,  stack rows below the main diagoal together,
      # the (i,j)-th element, for i>j, is in the position (i-2)*(i-1)/2+j
  # cat("STEP1","\n")

  ## Step 2: Apply FDR to identify r -- the number of connected pairs
  Ms=sort.int(M, index.return=T)
      # Ms$x are sorted P-values in ascending order, Ms$ix are the corresponding indices in M
  # cat(Ms$x, "\n\n")
  M=(Ms$x)*p0/beta
  # cat(M, "\n\n")
  Nn=1:p0
  r=Nn[M<=Nn]; j=length(r); if(j==0) stop("All component series are linearly independent")
  r=r[j]
  # cat("STEP2","\n")

  ## Step 3: find the pairs corresponding to the r maximum max_k |\rho(k)|
  h=mat.or.vec(p,1)
  for(i in 2:p) h[i]=(i-2)*(i-1)/2
  Inx=mat.or.vec(p,p); I=2:p
  for(k in 1:r) { 
    q=I[(Ms$ix[k]-h[I])>0]
    s=length(q)
    i=q[s]
    j=Ms$ix[k]-h[i]
    Inx[i,j]=1}
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
  if(K==0)
  output=list(NoGroups=0)
  else {
    Group=matrix(0,p,K)
    k=1
    for(j in 1:(p-1)) { if(Index[j]==1) { Group[,k]=G[,j]; k=k+1} }
    one_mem = which(!(c(1:p) %in%  Group))
    N2 = length(one_mem)
    if(N2>0)Group = cbind(Group,rbind(t(one_mem),matrix(0, p-1, N2)))[1:sum(Group[,1]>0),]
    else Group = as.matrix(Group[1:sum(Group[,1]>0),])
    q_block = K+N2
    Nosmem = c(N[N>0],rep(1,N2))
    if(verbose){
      cat("\n"); cat("Number of groups", q_block, "\n")
      cat("Number of members in those groups:", Nosmem, "\n")
      for(i in c(1:q_block)){
        cat("Groups",i,": contains these columns index of the zt:", drop(Group[,i]), "\n")
      }
    }
    output=list(NoGroups=q_block, No_of_Members=Nosmem, 
              Groups=Group)
  }
}

