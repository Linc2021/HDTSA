#' @name WN_FTS_test
#' @title Testing for functional white noise in high dimensions
#' @description Implements the high-dimensional white noise test for functional time series. 
#' It supports optional local linear smoothing and automatically switches between 
#' matrix and vectorized C++ implementations based on data dimensions.
#' 
#' @param Y_arr A \eqn{n \times p \times M} array, where \eqn{n} is the number of subjects, 
#' \eqn{p} is the dimension of the functional variables, and \eqn{M} is the number of time points.
#' @param U A \eqn{n \times M} matrix of observation time points. Required if \code{is.smooth = TRUE}.
#' @param lag.k The maximum time lag \eqn{L} to be tested. Default is 2.
#' @param B The number of bootstrap replications. Default is 1000.
#' @param alpha The significance level. Default is 0.05.
#' @param kernel.type The symmetric kernel used for long-run covariance estimation. 
#' Options: \code{"QS"} (default), \code{"Par"}, \code{"Bart"}.
#' @param is.smooth Logical. If \code{TRUE}, performs local linear smoothing on \code{Y_arr}. 
#' If \code{FALSE}, uses \code{Y_arr} directly. Default is \code{FALSE}.
#' 
#' @return A list of class \code{"hdtstest"}, which contains the following
#'   components:
#'
#'   \item{statistic}{The test statistic of the test.}
#'   \item{p.value}{The p-value of the test.}
#'   \item{crit.value}{The critical value of the test at the significance level alpha.} 
#'   \item{reject}{Reject or not.}
#'   \item{lag.k}{The lag for the test.}
#'   \item{kernel.type}{The type of kernel for the test.}
#'   \item{is.smooth}{Smooth the original data or not.}
#'   \item{computation.method}{Use matrix version to compute or vector version.}
#'      
#'   
#'   
#' @export
WN_FTS_test <- function(Y_arr, U = NULL, lag.k = 2, B = 1000, alpha = 0.05, 
                      kernel.type = c("QS", "Par", "Bart"),  is.smooth = FALSE) {
  
  # 1. Basic Parameter Preparation
  n <- dim(Y_arr)[1]
  p <- dim(Y_arr)[2]
  M <- dim(Y_arr)[3] # Number of time points
  
  kernel.type <- match.arg(kernel.type)
  ken_type <- switch(kernel.type, "QS" = 1, "Par" = 2, "Bart" = 3)
  
  # 2. Data Processing (Smoothing or Direct Use)
  if (is.smooth) {
    if (is.null(U)) stop("U (observation time points) must be provided if is.smooth = TRUE")
    
    # Target grid for smoothing (normalized to [0,1])
    uM <- matrix(rep(seq(0, 1, length = M), n), byrow = TRUE, nrow = n, ncol = M)
    
    # Perform local linear smoothing
    # local.linear.smooth handles CV for bandwidth and returns n*p*M array
    Vep_processed_arr <- local.linear.smooth(U, uM, Y_arr) 
  } else {
    # If data is already smooth or dense, use Y_arr directly
    Vep_processed_arr <- Y_arr
  }
  
  # 3. Reshape and Demean (Center the functional data)
  # Convert n*p*M array to n*(p*M) matrix
  Vep_temp <- t(apply(Vep_processed_arr, 1, function(x) as.vector(t(x))))
  # Global demeaning across subjects
  Vep <- Vep_temp - matrix(rep(colMeans(Vep_temp), n), nrow = n, byrow = TRUE)
  
  # 4. Computation Mode Selection
  # Switch to Vector version if p*M > 500 to save memory
  dim_size <- p * M
  use_vector <- dim_size > 500
  boot_nomal <- matrix(rnorm(B*(n-lag.k)), B, n-lag.k)
  
  if (use_vector) {
    # --- Vectorized Path (Memory Efficient) ---
    Tn <- TnC_vec(n, lag.k, Vep)
    bn <- WN_bandwith_vec(Vep, n, lag.k, p, M, ken_type)
    tt_boot <- newbootC_kernel_vec(n, lag.k, p, M, B, bn, ken_type, Vep, boot_nomal)
  } else {
    # --- Matrix Path (Speed Optimized for smaller dimensions) ---
    Tn <- TnC_mat(n, lag.k, Vep)
    # Matrix version requires explicit 'ft' matrix calculation
    ft <- WN_ftC_mat(n, lag.k, p, M, Vep) 
    bn <- WN_bandwith_mat(ft, n, lag.k, p, M, ken_type)
    tt_boot <- newbootC_kernel_mat(n, lag.k, p, M, B, bn, ken_type, Vep, boot_nomal)
  }
  
  # 5. Result Extraction
  tt_boot_sorted <- sort(tt_boot)
  p_value <- mean(tt_boot_sorted >= Tn)
  crit_value <- tt_boot_sorted[ceiling((1 - alpha) * B)]
  reject <- Tn > crit_value
  
  # 6. Format Output
  structure(list(
    statistic = Tn,
    p.value = p_value,
    crit.value = crit_value,
    reject = reject,
    lag.k = lag.k,
    kernel.type = kernel.type,
    is.smooth = is.smooth,
    computation.method = if(use_vector) "Vectorized" else "Matrix"
  ), class = "hdtstest")
}



#--------------------------------------1. Kernel smoothing-------------------------#
#------1.1 Kernel smoothing: kernel evaluation-----------------------#
# Function to calculate kernel weights for local linear smoothing
# U: Observed grid, uM: Target grid, bw: Bandwidth
kernel = function(U, uM, bw = 0.5){
  n = dim(uM)[1]     # Number of subjects
  M = dim(uM)[2]     # Number of target evaluation points
  Ntj = dim(U)[2]    # Number of observed time points
  
  # Ker_mat stores: [1] Kernel weight, [2] Weighted distance, [3] Weighted distance squared
  Ker_mat = array(0,dim = c(n,3,M,Ntj))
  
  for (i in 1:n){
    # Calculate the distance between each observed point and each target point
    U_delta = t(outer(U[i,],uM[i,],'-')) # M*Ntj matrix
    # Epanechnikov Kernel: 3/4 * (1 - u^2)
    Ker_mat[i,1,,] = 3 / (4 * bw) * (1 - (U_delta / bw)^2) * (abs(U_delta) < bw)
    # Basis for Local Linear Regression: (u-U_ij) and (u-U_ij)^2
    Ker_mat[i,2,,] = Ker_mat[i,1,,]*U_delta
    Ker_mat[i,3,,] = Ker_mat[i,1,,]*U_delta^2
  }
  
  return(Ker_mat)
}

#------1.2 Kernel smoothing: pre smoothing-----------------------#
# Function to estimate the smooth curve using local linear coefficients
smooth.curve.kernel = function(U,Y_arr,uM, Ker_mat){
  n = dim(uM)[1]
  M = dim(uM)[2]
  Ntj = dim(U)[2]
  p = dim(Y_arr)[2] # Number of variables 
  
  X = array(NA, dim = c(n,p,M))
  for (i in 1:n){
    # S matrix components: sum of K, K*delta, K*delta^2
    S = apply(Ker_mat[i,,,],c(1,2) ,sum)
    # Numerical stability check for the denominator of the local linear estimator
    denom_tmp = S[3,]*S[1,]-S[2,]^2
    denom = denom_tmp*(denom_tmp >= 10^(-4)) + 10^(-4)*(denom_tmp < 10^(-4))
    
    for (j in 1:p){
      # Weighted response calculation (T_0 and T_1)
      KZ = array(apply(Ker_mat[i,1:2,,], c(1),function(x) t(x)*Y_arr[i,j,] ),dim = c(Ntj,M,2))
      KZ = aperm(KZ, c(3, 2, 1))  # Reshape to 2*M*Ntj
      
      Tt = apply(KZ,c(1,2), sum) 
      # Local linear solution: intercept term represents the smoothed value at target point
      X[i,j,] =  (S[3,]*Tt[1,]-S[2,]*Tt[2,])/denom
    }
  }
  
  return(X)
} 

#------1.3 Cross-Validation for bandwidth -----------------------#
# Select optimal bandwidth using k-fold cross-validation
local.linear.smooth <- function(U, uM, Vep_arr, fold = 5){
  b_seq <- seq(0.1, 0.4, by=0.05) # Candidate bandwidths
  b_len <- length(b_seq)
  n <- dim(U)[1]
  Ntj <- dim(U)[2]
  M <- dim(uM)[2]
  p <- dim(Vep_arr)[2]
  
  residual <- matrix(NA, nrow=b_len, ncol=fold)
  flag_b <- 1
  for (bb in b_seq) {
    index <- sample(1:Ntj, replace = F)
    len <- ceiling(Ntj/fold)
    for (f in 1:fold) {
      # Splitting data into training and testing sets for CV
      if(f < fold){
        index_f <- index[((f-1)*len+1):(f*len)]
        sample_train <- Vep_arr[,,-index_f,drop=F]  
        sample_test <- Vep_arr[,,index_f,drop=F]
        U_train <- U[,-index_f,drop=F]
        U_test <- U[,index_f,drop=F]
        # Evaluate on test fold
        Ker_mat <- kernel(U_train, U_test, bb)
        sample_test_smooth <- smooth.curve.kernel(U=U_train, Y=sample_train, uM=U_test, Ker_mat=Ker_mat)
        residual[flag_b,f] <- mean((sample_test - sample_test_smooth)^2) 
      }else{
        index_f <- index[((f-1)*len+1):length(index)]
        sample_train <- Vep_arr[,,-index_f,drop=F]  
        sample_test <- Vep_arr[,,index_f,drop=F]
        U_train <- U[,-index_f,drop=F]
        U_test <- U[,index_f,drop=F]
        Ker_mat <- kernel(U_train, U_test, bb)
        sample_test_smooth <- smooth.curve.kernel(U=U_train, Y=sample_train, uM=U_test, Ker_mat=Ker_mat)
        residual[flag_b,f] <- mean((sample_test - sample_test_smooth)^2) 
      }
    }   # fold
    flag_b <- flag_b+1
  }  # b_seq
  
  # Pick bandwidth with minimum Mean Squared Error (MSE)
  b_opt <- b_seq[ which.min( apply(residual, 1, mean)) ]
  # Final smoothing with optimal bandwidth
  Ker_mat <- kernel(U, uM, b_opt)
  # Note: Uses full dataset 'U' and 'Vep_arr' for final estimation
  Vep_smooth_arr <- smooth.curve.kernel(U=U, Y_arr=Vep_arr, uM=uM, Ker_mat=Ker_mat) 
  return(Vep_smooth_arr)
}



