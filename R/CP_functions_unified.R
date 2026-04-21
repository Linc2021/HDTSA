#' @title Estimation of CP-factor models for high dimensional tensor-valued time series
#'
#' @description
#' \code{CP_MTS()} provides a unified interface for estimating CP-factor models
#' for tensor-valued time series. Let \eqn{\mathcal{Y}_t} be a tensor in
#' \eqn{\mathbb{R}^{p_1 \times \cdots \times p_m}}. The tensor CP-factor model is given by
#' \deqn{
#' \mathcal{Y}_t = \sum_{i=1}^d x_{t,i}\,\mathbf{a}_{i,1} \circ \mathbf{a}_{i,2} \circ \cdots \circ \mathbf{a}_{i,m} + \mathcal{E}_t, \quad t \ge 1,
#' }
#' where \eqn{1 \le d \le \min_{j \in [m]} p_j} is a fixed but unknown constant,
#' \eqn{\mathcal{E}_t \in \mathbb{R}^{p_1 \times \cdots \times p_m}} is the idiosyncratic error tensor,
#' \eqn{\mathbf{x}_t = (x_{t,1}, \ldots, x_{t,d})^\top} is the \eqn{d}-dimensional factor vector,
#' and \eqn{\mathbf{a}_{i,j}} is a \eqn{p_j}-dimensional loading vector
#' corresponding to the \eqn{i}-th factor and the \eqn{j}-th mode. Moreover, we assume \eqn{x_{t,i} = w_i f_{t,i}}
#' with the factor strength \eqn{w_i} and the standard factor process \eqn{f_{t,i}}.
#'
#' When \eqn{m=2}, the model reduces to the matrix time series CP-factor model discussed
#' in Chang et al. (2023):
#' \deqn{
#' \mathbf{Y}_t = \mathbf{A}\mathbf{X}_t\mathbf{B}^\top + \boldsymbol{\epsilon}_t,
#' }
#' where \eqn{\mathbf{Y}_t \in \mathbb{R}^{p \times q}} with \eqn{p = p_1} and \eqn{q = p_2},
#' \eqn{\mathbf{X}_t = \mathrm{diag}(x_{t,1}, \ldots, x_{t,d})},
#' \eqn{\boldsymbol{\epsilon}_t \in \mathbb{R}^{p \times q}},
#' \eqn{\mathbf{A} = (\mathbf{a}_{1,1}, \ldots, \mathbf{a}_{d,1}) \in \mathbb{R}^{p \times d}},
#' and \eqn{\mathbf{B} = (\mathbf{a}_{1,2}, \ldots, \mathbf{a}_{d,2}) \in \mathbb{R}^{q \times d}}
#' are loading matrices.
#'
#' This function supports several estimation methods for tensor CP-factor models.
#' When the tensor time series is of order \eqn{m>2}, the available method is
#' \code{"CP.DPI"}. When \eqn{m=2}, that is, when the tensor-valued time series
#' reduces to a matrix-valued time series, the available methods are
#' \code{"CP.Direct"}, \code{"CP.Refined"}, \code{"CP.Unified"}, and
#' \code{"CP.DPI"}. The first three follow Chang et al. (2023, 2026a+), while
#' \code{"CP.DPI"} is implemented through Chang et al. (2026b+).
#'
#' @details
#' The function is designed for tensor-valued time series stored as an array of
#' dimension \code{c(n, p1, ..., pm)}. When \code{dim(Y)} has length greater than 3,
#' the observed process is a tensor-valued time series of order \eqn{m>2}, and
#' only \code{method = "CP.DPI"} is available. When \code{dim(Y)} has length 3,
#' the observed process is a second-order tensor time series, that is, a
#' matrix-valued time series corresponding to the special case \eqn{m=2}. In this
#' case, all four methods, namely \code{"CP.Direct"}, \code{"CP.Refined"},
#' \code{"CP.Unified"}, and \code{"CP.DPI"}, can be used.
#'
#' When the tensor-valued time series is of order \eqn{m=2}, the model reduces to
#' the matrix CP-factor model, and the methods \code{"CP.Direct"},
#' \code{"CP.Refined"}, and \code{"CP.Unified"} are available. All three methods
#' involve the estimation of the autocovariance between \eqn{\mathbf{Y}_t} and
#' \eqn{\xi_t} at lag \eqn{k}, which is defined as follows:
#' \deqn{
#' \hat{\mathbf{\Sigma}}_{k}
#' = T_{\delta_1}\{\hat{\boldsymbol{\Sigma}}_{\mathbf{Y},\xi}(k)\}
#' \quad \mbox{with} \quad
#' \hat{\boldsymbol{\Sigma}}_{\mathbf{Y},\xi}(k)
#' = \frac{1}{n-k}\sum_{t=k+1}^n(\mathbf{Y}_t-\bar{\mathbf{Y}})(\xi_{t-k}-\bar{\xi}),
#' }
#' where \eqn{\bar{\mathbf{Y}} = n^{-1}\sum_{t=1}^n \mathbf{Y}_t},
#' \eqn{\bar{\xi} = n^{-1}\sum_{t=1}^n \xi_t}, and \eqn{T_{\delta_1}(\cdot)}
#' is a threshold operator defined as
#' \eqn{T_{\delta_1}(\mathbf{W}) = \{w_{i,j}1(|w_{i,j}| \ge \delta_1)\}}
#' for any matrix \eqn{\mathbf{W}=(w_{i,j})}, with threshold level
#' \eqn{\delta_1 \ge 0} and \eqn{1(\cdot)} denoting the indicator function.
#' Chang et al. (2024) and Chang et al. (2026a+) suggest choosing
#' \eqn{\delta_1 = 0} when \eqn{p_1} and \eqn{p_2} are fixed, and
#' \eqn{\delta_1 > 0} when \eqn{p_1p_2 \gg n}.
#'
#' The refined estimation method further involves
#' \deqn{
#' \check{\mathbf{\Sigma}}_{k}
#' = T_{\delta_2}\{\hat{\mathbf{\Sigma}}_{\check{\mathbf{Y}}}(k)\}
#' \quad \mbox{with} \quad
#' \hat{\mathbf{\Sigma}}_{\check{\mathbf{Y}}}(k)
#' = \frac{1}{n-k}\sum_{t=k+1}^n(\mathbf{Y}_t-\bar{\mathbf{Y}})
#' \otimes \mathrm{vec}(\mathbf{Y}_{t-k}-\bar{\mathbf{Y}}),
#' }
#' where \eqn{T_{\delta_2}(\cdot)} is a threshold operator with threshold level
#' \eqn{\delta_2 \ge 0}, and \eqn{\mathrm{vec}(\cdot)} is the vectorization operator,
#' with \eqn{\mathrm{vec}(\mathbf{H})} being the \eqn{(m_1m_2)\times 1} vector
#' obtained by stacking the columns of the \eqn{m_1 \times m_2} matrix \eqn{\mathbf{H}}.
#' See Section 3.2.2 of Chang et al. (2024) for details.
#'
#' The unified estimation method further involves
#' \deqn{
#' \vec{\mathbf{\Sigma}}_{k}
#' = T_{\delta_3}\{\hat{\boldsymbol{\Sigma}}_{\mathrm{vec}(\mathbf{Y})}(k)\}
#' \quad \mbox{with} \quad
#' \hat{\boldsymbol{\Sigma}}_{\mathrm{vec}(\mathbf{Y})}(k)
#' = \frac{1}{n-k}\sum_{t=k+1}^n
#' \mathrm{vec}(\mathbf{Y}_t-\bar{\mathbf{Y}})
#' \{\mathrm{vec}(\mathbf{Y}_{t-k}-\bar{\mathbf{Y}})\}',
#' }
#' where \eqn{T_{\delta_3}(\cdot)} is a threshold operator with threshold level
#' \eqn{\delta_3 \ge 0}. See Section 4.2 of Chang et al. (2026a+) for details.
#'
#' For tensor-valued time series of general order \eqn{m \ge 2}, the method
#' \code{"CP.DPI"} is also available. This method is based on the Double Projection
#' Iterations (DPI) procedure for the tensor CP-factor model; see Chang et al. (2026b+)
#' for details. In particular, when \eqn{m>2}, \code{"CP.DPI"} is the only
#' available method, while when \eqn{m=2}, it can be used together with
#' \code{"CP.Direct"}, \code{"CP.Refined"}, and \code{"CP.Unified"}.
#'
#' @param Y An array representing a tensor-valued time series, with dimension
#'   \code{c(n, p1, ..., pm)}, where \eqn{n} is the sample size. In particular,
#'   when \code{dim(Y)} has length 3, \code{Y} corresponds to a matrix-valued
#'   time series, which is the special case \eqn{m=2}.
#'
#' @param xi An auxiliary scalar series \eqn{(\xi_1,\ldots,\xi_n)^\top}, which is
#'   a linear combination of \eqn{\mathrm{vec}(\mathcal{Y}_t)}. If \code{xi = NULL}
#'   (the default) and \code{method} is one of \code{"CP.Direct"},
#'   \code{"CP.Refined"}, or \code{"CP.Unified"}, \eqn{\xi_t} is determined by
#'   the PCA method introduced in Section 5.1 of Chang et al. (2024). If
#'   \code{xi = NULL} and \code{method = "CP.DPI"}, then \eqn{\xi_t} can be
#'   estimated by the PCA method or by a random projection method by setting
#'   \code{Random.Projection = TRUE} in \code{control.DPI}.
#'
#' @param Rank The prescribed number of factors. \code{Rank} is a list containing
#'   the relevant structural ranks, such as \code{d}, \code{d1}, and \code{d2},
#'   depending on the method. If set to \code{NULL} (the default), these quantities
#'   are estimated from the data. If \code{method = "CP.DPI"}, the number of
#'   factors can be selected by the ER method or the Log-ER method by setting
#'   \code{Ratio.type} in \code{control.DPI} to \code{"classical"} or \code{"log"},
#'   respectively.
#'
#' @param lag.k The lag parameter \eqn{K} used in \code{"CP.Refined"} and
#'   \code{"CP.Unified"}. When \eqn{m=2}, it is used to construct the
#'   nonnegative definite matrices \eqn{\hat{\mathbf{M}}_1} and
#'   \eqn{\hat{\mathbf{M}}_2}:
#'   \deqn{
#'   \hat{\mathbf{M}}_1 = \sum_{k=1}^{K}\hat{\mathbf{\Sigma}}_{k}\hat{\mathbf{\Sigma}}_{k}'
#'   \quad \mbox{and} \quad
#'   \hat{\mathbf{M}}_2 = \sum_{k=1}^{K}\hat{\mathbf{\Sigma}}_{k}'\hat{\mathbf{\Sigma}}_{k},
#'   }
#'   where \eqn{\hat{\mathbf{\Sigma}}_{k}} is the estimated cross-covariance
#'   between \eqn{\mathbf{Y}_t} and \eqn{\xi_t} at lag \eqn{k}. The default is 20.
#'
#' @param lag.ktilde The lag parameter \eqn{\tilde K} involved in the unified
#'   estimation method [see (16) in Chang et al. (2026a+)], which is used only when
#'   \code{method = "CP.Unified"}. The default is 10.
#'
#' @param method A character string specifying the estimation method. Available
#'   choices are \code{"CP.Direct"}, \code{"CP.Refined"}, \code{"CP.Unified"},
#'   and \code{"CP.DPI"}. When \eqn{m=2}, all four methods can be used. When
#'   \eqn{m>2}, only \code{"CP.DPI"} is available. The validity of
#'   \code{"CP.Direct"} and \code{"CP.Refined"} depends on the assumption
#'   \eqn{d_1=d_2=d}. When \eqn{d_1,d_2 \le d}, \code{"CP.Unified"} can be applied.
#'   See Chang et al. (2024, 2026a+) for details.
#'
#' @param thresh1 Logical. If \code{FALSE} (the default), no thresholding is
#'   applied in \eqn{\hat{\mathbf{\Sigma}}_{k}}, corresponding to \eqn{\delta_1=0}.
#'   If \code{TRUE}, \eqn{\delta_1} is set through \code{delta1}. This argument is
#'   used by \code{"CP.Direct"}, \code{"CP.Refined"}, and \code{"CP.Unified"}.
#'
#' @param thresh2 Logical. If \code{FALSE} (the default), no thresholding is
#'   applied in \eqn{\check{\mathbf{\Sigma}}_{k}}, corresponding to \eqn{\delta_2=0}.
#'   If \code{TRUE}, \eqn{\delta_2} is set through \code{delta2}. This argument is
#'   used only when \code{method = "CP.Refined"}.
#'
#' @param thresh3 Logical. If \code{FALSE} (the default), no thresholding is
#'   applied in \eqn{\vec{\mathbf{\Sigma}}_{k}}, corresponding to \eqn{\delta_3=0}.
#'   If \code{TRUE}, \eqn{\delta_3} is set through \code{delta3}. This argument is
#'   used only when \code{method = "CP.Unified"}.
#'
#' @param delta1 The threshold level \eqn{\delta_1}, used by
#'   \code{"CP.Direct"}, \code{"CP.Refined"}, and \code{"CP.Unified"}.
#'   The default is \eqn{2\sqrt{n^{-1}\log(p_1p_2)}}.
#'
#' @param delta2 The threshold level \eqn{\delta_2}, used only when
#'   \code{method = "CP.Refined"}. The default is
#'   \eqn{2\sqrt{n^{-1}\log(p_1p_2)}}.
#'
#' @param delta3 The threshold level \eqn{\delta_3}, used only when
#'   \code{method = "CP.Unified"}. The default is
#'   \eqn{2\sqrt{n^{-1}\log(p_1p_2)}}.
#'
#' @param A.init Optional initial loading matrices used only when
#'   \code{method = "CP.DPI"}. It should be a list of length
#'   \eqn{m}, where the \eqn{j}-th element is a
#'   \eqn{p_j \times d} matrix. If \code{NULL}, an initial estimator is obtained
#'   from an initialization step.
#'
#' @param control.DPI A named list of control parameters used when
#'   \code{method = "CP.DPI"}. The supported components are:
#'   \describe{
#'     \item{\code{lag.k.DPI}}{Positive integer. Number of lags used in
#'       \eqn{\tilde{\mathbf{M}}_j}. Default is \eqn{10}.}
#'     \item{\code{Threshold}}{Logical. Whether thresholding is applied in the
#'       initialization step. Default is \code{TRUE}.}
#'     \item{\code{delta}}{Optional thresholding level used in the initialization
#'       step. Default is \code{NULL}. If \code{NULL}, it is selected via a grid
#'       search method.}
#'     \item{\code{delta2}}{Numeric vector of length \code{m}, controlling the
#'       thresholding level in each tensor mode during the iterative update. The
#'       default \code{j}-th element is
#'       \eqn{\hat{\sigma}_0 (n^{-1}\log p_j)^{1/2}}, where
#'       \eqn{\hat{\sigma}_0^2 = (n \prod_{j=1}^m p_j)^{-1}
#'       \sum_{t=1}^n\|\mathcal{Y}_t\|^2_{\mathrm{F}}}.}
#'     \item{\code{Ratio.type}}{Character string specifying the ratio criterion
#'       used in rank estimation. Typical choices are \code{"log"} (Log-ER method)
#'       and \code{"classical"} (ER method). Default is \code{"log"}.}
#'     \item{\code{Random.Projection}}{Logical. If \code{TRUE}, a randomized
#'       projection step is used to select \code{xi}. Default is \code{FALSE}.}
#'     \item{\code{iter_max}}{Maximum number of iterative updates. Default is
#'       \eqn{20}.}
#'     \item{\code{eps}}{Stopping tolerance for the iterative algorithm.
#'       Default is \eqn{10^{-4}}.}
#'     \item{\code{grid.num}}{Integer. Number of grid points used when selecting
#'       the thresholding level in the initialization step. Default is \eqn{50}.}
#'     \item{\code{delta_max}}{Maximum value of the thresholding grid in the
#'       initialization step. Default is \eqn{0.1 \hat{\sigma}_0 (n^{-1} \sum_{j = 1}^m \log p_j)^{1/2}}.}
#'     \item{\code{print.eps}}{Logical. Whether to print the iterative
#'       convergence measure. Default is \code{FALSE}.}
#'     \item{\code{iter_lag}}{Positive integer. Number of candidate lags used in
#'       each iterative update. Default is \eqn{1}.}
#'     \item{\code{all.put}}{Logical. If \code{TRUE}, the iterative routine
#'       returns full intermediate outputs; otherwise only a compact result is
#'       returned. Default is \code{FALSE}.}
#'     \item{\code{A}}{Optional true loading matrices, used only for diagnostic
#'       purposes in simulations. Default is \code{NULL}.}
#'     \item{\code{Component}}{Optional true common component tensor, used only
#'       for diagnostic purposes such as CP reconstruction loss in simulations.
#'       Default is \code{NULL}.}
#'   }
#'
#' @return
#' When \code{method = "CP.DPI"}, the function returns a list with three components:
#' \describe{
#'   \item{\code{res.iter}}{A list containing the iterative loading matrices
#'   \code{A.hat}, the initial loading matrices \code{A.init}, the estimated factor
#'   series \code{f_hat}, the initial factor series \code{f_hat_inl}, the number
#'   of iterations \code{iter_step}, and the thresholded moment matrices used in
#'   inference \code{Sigma.yij.xii.1}.}
#'   \item{\code{res.init}}{The initial estimator when \code{A.init = NULL}. If
#'   \code{A.init} is supplied by the user, this component is \code{NULL}.}
#'   \item{\code{control.DPI}}{The control list actually used in the function
#'   after merging user-supplied values with the defaults.}
#' }
#'
#' When \eqn{m=2} and one of \code{"CP.Direct"}, \code{"CP.Refined"}, or
#' \code{"CP.Unified"} is used, the function returns an object of class
#' \code{"mtscp"} containing the following components:
#' \describe{
#'   \item{\code{A}}{The estimated \eqn{p_1 \times \hat{d}} loading matrix
#'   corresponding to the first mode.}
#'   \item{\code{B}}{The estimated \eqn{p_2 \times \hat{d}} loading matrix
#'   corresponding to the second mode.}
#'   \item{\code{f}}{The estimated latent process
#'   \eqn{\hat{\mathbf{x}}_t=(\hat{x}_{t,1},\ldots,\hat{x}_{t,\hat{d}})^\top}.}
#'   \item{\code{Rank}}{The estimated rank information.}
#'   \item{\code{method}}{A string indicating which estimation method is used.}
#' }
#'
#' @examples
#' # Example 1: matrix-valued time series with d = d1 = d2.
#' p <- 10
#' q <- 10
#' n <- 400
#' d = d1 = d2 <- 3
#' ## DGP.CP() generates simulated data for the example in Chang et al. (2026a+).
#' data <- DGP.CP(n, p, q, d, d1, d2)
#' Y <- data$Y
#'
#' ## d is unknown
#' res1 <- CP_MTS(Y, method = "CP.Direct")
#' res2 <- CP_MTS(Y, method = "CP.Refined")
#' res3 <- CP_MTS(Y, method = "CP.Unified")
#' res4 <- CP_MTS(Y, method = "CP.DPI")
#'
#' ## d is known
#' res5 <- CP_MTS(Y, Rank = list(d = 3), method = "CP.Direct")
#' res6 <- CP_MTS(Y, Rank = list(d = 3), method = "CP.Refined")
#' res7 <- CP_MTS(Y, Rank = list(d = 3), method = "CP.DPI")
#'
#'
#' # Example 2: matrix-valued time series with d1, d2 < d.
#' p <- 10
#' q <- 10
#' n <- 400
#' d1 = d2 <- 2
#' d <- 3
#' data <- DGP.CP(n, p, q, d, d1, d2)
#' Y1 <- data$Y
#'
#' ## d, d1 and d2 are unknown
#' res8 <- CP_MTS(Y1, method = "CP.Unified")
#'
#' ## d, d1 and d2 are known
#' res9 <- CP_MTS(Y1, Rank = list(d = 3, d1 = 2, d2 = 2), method = "CP.Unified")
#'
#'
#'
#' @references
#' Chang, J., He, J., Yang, L., & Yao, Q. (2023). Modelling matrix time series
#' via a tensor CP-decomposition. \emph{Journal of the Royal Statistical Society
#' Series B: Statistical Methodology}, \strong{85}, 127--148.
#' \url{doi:10.1093/jrsssb/qkac011}.
#'
#' Chang, J., Du, Y., Huang, G., & Yao, Q. (2026a+). Identification and
#' estimation for matrix time series CP-factor models. \emph{The Annals of
#' Statistics}, in press. \url{doi:10.48550/arXiv.2410.05634}.
#'
#' Chang, J., Huang, G., Yao, Q., & Yu, L. (2026b+). CP-Factorization for High Dimensional Tensor Time
#' Series and Double Projection Iterations. \emph{Journal of the Royal Statistical Society
#' Series B: Statistical Methodology}, major revision.
#'
#' @export
#' @useDynLib HDTSA
#' @importFrom stats arima.sim rnorm runif rt

CP_MTS = function(Y,
                  xi = NULL,
                  Rank = NULL,
                  lag.k = 20,
                  lag.ktilde = 10,
                  method = c("CP.Direct", "CP.Refined", "CP.Unified", "CP.DPI"),
                  thresh1 = FALSE,
                  thresh2 = FALSE,
                  thresh3 = FALSE,
                  delta1 = 2 * sqrt(log(dim(Y)[2] * dim(Y)[3]) / dim(Y)[1]),
                  delta2 = delta1,
                  delta3 = delta1,
                  A.init = NULL,
                  control.DPI = list()
) {
  
  method <- match.arg(method)
  ydim <- dim(Y)
  if (is.null(ydim) || length(ydim) < 3) {
    stop("Y must be an array of dimension c(n, p1, ..., pm) with m >= 2.")
  }
  
  m <- length(ydim) - 1
  
  ## =========================================================
  ## Tensor CP-factor model: DPI method
  ## =========================================================
  if (method == "CP.DPI") {
    if(!is.null(Rank)){
      Rank = Rank$d
    }
    
    return(
      HDTTS.CP.iter.DPI(
        Y = Y,
        A.init = A.init,
        xi = xi,
        Rank = Rank,
        control.DPI = control.DPI
      )
    )
  }
  
  ## =========================================================
  ## Matrix CP-factor model methods
  ## =========================================================
  if (m != 2) {
    stop("Methods 'CP.Direct', 'CP.Refined', and 'CP.Unified' are only available when m = 2. For m > 2, please use method = 'CP.DPI'.")
  }
  
  n <- ydim[1]
  p <- ydim[2]
  q <- ydim[3]
  
  if (is.null(xi)) {
    xi <- est.xi(Y)$xi
  }
  
  if(method == "CP.Direct"){
    S_yxi_1 <- Autocov_xi_Y(Y, xi, k = 1, thresh = thresh1, delta = delta1)
    S_yxi_2 <- Autocov_xi_Y(Y, xi, k = 2, thresh = thresh1, delta = delta1)
    
    if(p > q){
      K1 <- t(S_yxi_1) %*% S_yxi_1
      eg1 <- eigen(K1)
      w <- eg1$values
      ww <- w[-1] / w[-length(w)]
      d <- which(ww == min(ww[1:floor(0.75 * q)]))
      
      if(!is.null(Rank)){
        if(!is.null(Rank$d)){
          d <- Rank$d
        } else {
          stop("List Rank without d, use Rank = list(d = ?)")
        }
      }
      
      if (d > 1){
        K1 <- eg1$vectors[, 1:d] %*% diag(eg1$values[1:d]) %*% t(eg1$vectors[, 1:d])
      } else {
        K1 <- eg1$vectors[, 1] %*% diag(eg1$values[1], 1) %*% t(eg1$vectors[, 1])
      }
      K2 <- t(S_yxi_1) %*% S_yxi_2
      
      Geg <- geigen::geigen(K2, K1)
      evalues <- Geg$values[which(Mod(Geg$values) <= 10^5 & Geg$values != 0)]
      
      Bl <- Geg$vectors[, which(Geg$values %in% evalues), drop = FALSE]
      A <- apply(S_yxi_1 %*% Bl, 2, l2s)
      Al <- t(MASS::ginv(A))
      B <- apply(t(S_yxi_1) %*% Al, 2, l2s)
    } else {
      K1 <- S_yxi_1 %*% t(S_yxi_1)
      eg1 <- eigen(K1)
      w <- eg1$values
      ww <- w[-1] / w[-length(w)]
      d <- which(ww == min(ww[1:floor(0.75 * p)]))
      
      if(!is.null(Rank)){
        if(!is.null(Rank$d)){
          d <- Rank$d
        } else {
          stop("List Rank without d, use Rank = list(d = ?)")
        }
      }
      
      if (d > 1){
        K1 <- eg1$vectors[, 1:d] %*% diag(eg1$values[1:d]) %*% t(eg1$vectors[, 1:d])
      } else {
        K1 <- eg1$vectors[, 1] %*% diag(eg1$values[1], 1) %*% t(eg1$vectors[, 1])
      }
      K2 <- S_yxi_1 %*% t(S_yxi_2)
      
      Geg <- geigen::geigen(K2, K1)
      evalues <- Geg$values[which(Mod(Geg$values) <= 10^5 & Geg$values != 0)]
      Al <- Geg$vectors[, which(Geg$values %in% evalues), drop = FALSE]
      B <- apply(t(S_yxi_1) %*% Al, 2, l2s)
      Bl <- t(MASS::ginv(B))
      A <- apply(S_yxi_1 %*% Bl, 2, l2s)
    }
    
    H <- matrix(NA, p * q, d)
    for (ii in 1:d) {
      H[, ii] <- B[, ii] %x% A[, ii]
    }
    f <- Vec.tensor(Y) %*% H %*% MASS::ginv(t(H) %*% H)
    
    if(is.complex(A) || is.complex(B)){
      A <- Complex2Real(A)
      B <- Complex2Real(B)
      f <- Complex2Real(f)
    }
    
    con <- structure(
      list(A = A, B = B, f = f, Rank = list(d = d), method = method),
      class = "mtscp"
    )
    return(con)
  }
  
  if(method == "CP.Refined"){
    M1 = M2 <- 0
    dmax <- round(min(p, q) * 0.75)
    
    for (kk in 1:lag.k){
      S_yxi_k <- Autocov_xi_Y(Y, xi, k = kk, thresh = thresh1, delta = delta1)
      M1 <- M1 + S_yxi_k %*% t(S_yxi_k)
      M2 <- M2 + t(S_yxi_k) %*% S_yxi_k
    }
    
    ev_M1 <- eigen(M1)
    ev_M2 <- eigen(M2)
    
    d1 <- which.max(ev_M1$values[1:dmax] / ev_M1$values[2:(dmax + 1)])
    d2 <- which.max(ev_M2$values[1:dmax] / ev_M2$values[2:(dmax + 1)])
    d <- ifelse(p > q, d1, d2)
    
    if(!is.null(Rank)){
      if(!is.null(Rank$d)){
        d <- Rank$d
      } else {
        stop("List Rank without d, use Rank = list(d = ?)")
      }
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
    } else {
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
      } else {
        S_Zxi_1 <- Autocov_xi_Y(Y = Z, xi = xi$xi, k = 1)
        S_Zxi_2 <- Autocov_xi_Y(Y = Z, xi = xi$xi, k = 2)
      }
      
      vl <- eigen(MASS::ginv(t(S_Zxi_1) %*% S_Zxi_1) %*% t(S_Zxi_1) %*% S_Zxi_2)$vectors
      ul <- eigen(MASS::ginv(S_Zxi_1 %*% t(S_Zxi_1)) %*% S_Zxi_1 %*% t(S_Zxi_2))$vectors
      
      U <- apply(S_Zxi_1 %*% vl, 2, l2s)
      V <- apply(t(S_Zxi_1) %*% ul, 2, l2s)
      
      A <- P %*% U
      B <- Q %*% V
      
      W <- matrix(NA, d^2, d)
      for (ii in 1:d) {
        W[, ii] <- V[, ii] %x% U[, ii]
      }
      f <- Vec.tensor(Z) %*% W %*% solve(t(W) %*% W)
      
      if(is.complex(A) || is.complex(B)){
        A <- Complex2Real(A)
        B <- Complex2Real(B)
        f <- Complex2Real(f)
      }
    }
    
    con <- structure(
      list(A = A, B = B, f = f, Rank = list(d = d), method = method),
      class = "mtscp"
    )
    return(con)
  }
  
  if(method == "CP.Unified"){
    if(is.null(Rank)){
      PQ_hat_tol <- est.d1d2.PQ(Y, xi, K = lag.k, thresh = thresh1, delta = delta1)
      d1 <- PQ_hat_tol$d1_hat
      d2 <- PQ_hat_tol$d2_hat
      d  <- NULL
      P  <- PQ_hat_tol$P_hat
      Q  <- PQ_hat_tol$Q_hat
      
      if(d1 == 1 || d2 == 1) d <- d1 * d2
    } else {
      if(all(!is.null(Rank$d), !is.null(Rank$d1), !is.null(Rank$d2))){
        d  <- Rank$d
        d1 <- Rank$d1
        d2 <- Rank$d2
      } else {
        stop("List Rank without d, d1 and d2, use Rank = list(d = ?, d1 = ?, d2 = ?)")
      }
      
      PQ_hat_tol <- est.PQ(Y, xi, d1, d2, K = lag.k, thresh = thresh1, delta = delta1)
      P <- PQ_hat_tol$P_hat
      Q <- PQ_hat_tol$Q_hat
    }
    
    if(d1 == 1 & d2 == 1){
      d <- 1
      f <- vector()
      for (tt in 1:n) {
        f[tt] <- t(P) %*% Y[tt, , ] %*% Q
      }
      rank <- list(d = d, d1 = d1, d2 = d2)
      con <- structure(
        list(A = as.matrix(P), B = as.matrix(Q), f = as.matrix(f), Rank = rank, method = method),
        class = "mtscp"
      )
      return(con)
    } else {
      if(is.null(d)){
        W_hat_tol <- est.d.Wf(Y, P, Q, Ktilde = lag.ktilde, thresh = thresh3, delta = delta3)
        d <- W_hat_tol$d_hat
        W <- W_hat_tol$W_hat
        f <- W_hat_tol$f_hat
      } else {
        W_hat_tol <- est.Wf(Y, P, Q, d, Ktilde = lag.ktilde, thresh = thresh3, delta = delta3)
        W <- W_hat_tol$W_hat
        f <- W_hat_tol$f_hat
      }
      
      if(d1 == 1 || d2 == 1){
        Theta <- NULL
        if(d1 == 1){
          U <- 1
          V <- W
          stop("d1 equal to 1, V cannot be identified uniquely!")
        }
        if(d2 == 1){
          U <- W
          V <- 1
          stop("d2 equal to 1, U cannot be identified uniquely!")
        }
        if(d1 == 1 & d2 == 1){
          U <- 1
          V <- 1
        }
        U <- as.matrix(U)
        V <- as.matrix(V)
      } else {
        UV_hat_tol <- est.UV.JAD(W, d1, d2, d)
        U <- UV_hat_tol$U
        V <- UV_hat_tol$V
        Theta <- UV_hat_tol$Theta
      }
      
      A <- P %*% U
      B <- Q %*% V
      rank <- list(d = d, d1 = d1, d2 = d2)
      con <- structure(
        list(A = as.matrix(A), B = as.matrix(B), f = as.matrix(f), Rank = rank, method = method),
        class = "mtscp"
      )
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
  
  real.index <- which(colSums(abs(IMA)) == 0)
  
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


HDTTS.CP.iter.DPI <- function(Y,
                              A.init = NULL,
                              xi = NULL,
                              Rank = NULL,
                              control.DPI = list()
) {
  n <- dim(Y)[1]
  D <- dim(Y)[-1]
  m <- length(D)
  
  control.DPI.default <- list(
    lag.k.DPI = 10,
    Threshold = TRUE,
    delta = NULL,
    delta2 = rep(1, length(dim(Y)) - 1),
    Ratio.type = "log",
    Random.Projection = FALSE,
    iter_max = 20,
    eps = 1e-4,
    grid.num = 50,
    delta_max = 0.1,
    print.eps = FALSE,
    iter_lag = 1,
    all.put = FALSE,
    A = NULL,
    Component = NULL
  )
  
  control.DPI <- utils::modifyList(control.DPI.default, control.DPI)
  
  if (!is.null(control.DPI$lag.k) && is.null(control.DPI$lag.k.DPI)) {
    control.DPI$lag.k.DPI <- control.DPI$lag.k
  }
  if (!is.null(control.DPI$Random.Project) && is.null(control.DPI$Random.Projection)) {
    control.DPI$Random.Projection <- control.DPI$Random.Project
  }
  
  lag.k <- control.DPI$lag.k.DPI
  Threshold <- control.DPI$Threshold
  delta <- control.DPI$delta
  delta2 <- control.DPI$delta2
  Ratio.type <- control.DPI$Ratio.type
  Random.Projection <- control.DPI$Random.Projection
  iter_max <- control.DPI$iter_max
  eps <- control.DPI$eps
  grid.num <- control.DPI$grid.num
  delta_max <- control.DPI$delta_max
  print.eps <- control.DPI$print.eps
  iter_lag <- control.DPI$iter_lag
  all.put <- control.DPI$all.put
  A <- control.DPI$A
  Component <- control.DPI$Component
  
  if (length(delta2) != m) {
    stop("length(delta2) must equal length(dim(Y)) - 1.")
  }
  
  if (is.null(xi)) {
    xi.chang <- tensor.est.xi(Y)
    
    if (isTRUE(Random.Projection)) {
      res.only.used.rank <- HDTTS.CP.est(
        Y = Y,
        xi = xi.chang,
        K = lag.k,
        Ratio.type = Ratio.type,
        grid_delta1 = grid.num,
        delta_max = delta_max
      )
      
      r_hat_first <- res.only.used.rank$r.hat
      r_breve <- 2 * r_hat_first
      if (!is.null(Rank)) {
        r_breve <- 2 * Rank
      }
      
      xi.res <- RP.xi.sel(
        Y = Y,
        r_breve = r_breve,
        eps = 0.1,
        lag.k = lag.k,
        Randomized.time = 50,
        A = A
      )
      xi <- xi.res$xi.sel
    } else {
      xi <- xi.chang
    }
  }
  
  if (is.null(A.init)) {
    res.init <- HDTTS.CP.est(
      Y = Y,
      xi = xi,
      Rank = Rank,
      K = lag.k,
      Threshold = Threshold,
      delta = delta,
      Ratio.type = Ratio.type,
      grid_delta1 = grid.num,
      delta_max = delta_max
    )
    
    A.hat <- res.init$A.hat
    sigma0 <- res.init$sigma0
  } else {
    res.init <- NULL
    A.hat <- A.init
    sigma0 <- sqrt(sum(Y^2) / (n * prod(D)))
  }
  
  if (isFALSE(Threshold)) {
    delta2 <- rep(0, m)
  }
  
  res.iter.fin <- CP.iter.DPI.xi(
    Y = Y,
    A.hat = A.hat,
    K = lag.k,
    n = n,
    delta2 = delta2,
    sigma0 = sigma0,
    iter_max = iter_max,
    eps = eps,
    print.eps = print.eps,
    iter_lag = iter_lag,
    all.put = all.put,
    A = A,
    Component = Component
  )
  
  list(
    res.iter = res.iter.fin,
    res.init = res.init,
    control.DPI = control.DPI
  )
}


#' Inference for the estimates obtained by Double Projection Iterations (DPI) method in the tensor CP-factor model
#'
#' This function performs inference for a linear functional of the debiased
#' DPI estimator of the factor loading vector in the tensor CP-factor model.
#' Given a direction vector \code{h}, a factor index \code{i}, a mode index \code{j},
#' the observed tensor series \code{Y}, and the output object \code{res.CP.DPI} of the DPI method,
#' the function returns the debiased estimate of \eqn{\mathbf{h}^\top (\hat{\mathbf{a}}_{i,j} - \hat{\mathbf{\vartheta}}_{i.j})},
#' its estimated standard error, the corresponding DPI estimator \eqn{\hat{\mathbf{a}}_{i,j}}, and the
#' estimated bias-correction term \eqn{\hat{\mathbf{\vartheta}}_{i.j})}.
#'
#' @param h A numeric vector of length \eqn{d_j}. This specifies the linear
#'   functional \eqn{\mathbf{h}^\top (\hat{\mathbf{a}}_{i,j} - \hat{\mathbf{\vartheta}}_{i.j})} of interest.
#' @param i A positive integer. The factor index.
#' @param j A positive integer. The tensor mode index.
#' @param Y An array containing the observed tensor time series. Its dimension is
#'   typically \eqn{n \times d_1 \times \cdots \times d_m}.
#' @param res.CP.DPI An output object from \code{CP_MTS} with \code{method = "CP.DPI"}.
#'
#'
#' @details
#' The reported standard error is
#' \deqn{
#' \sqrt{\widehat{\mathrm{Var}} \big\{\mathbf{h}^\top (\hat{\mathbf{a}}_{i,j} - \hat{\mathbf{\vartheta}}_{i.j})\big\}/n} \,,
#' }
#' where \eqn{n} is the sample size.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{aij.h.de}}{ \eqn{\mathbf{h}^\top (\hat{\mathbf{a}}_{i,j} - \hat{\mathbf{\vartheta}}_{i.j})}.}
#'   \item{\code{se.h.ij}}{The estimated standard error of \code{aij.h.de}.}
#'   \item{\code{aij.iter}}{The original DPI estimator of the loading vector \eqn{\mathbf{a}_{i,j}}.}
#'   \item{\code{vartheta.ij}}{The estimated bias-correction term \eqn{\hat{\mathbf{\vartheta}}_{i.j}} used in debiasing.}
#' }
#'
#'
#' @examples
#' \dontrun{
#' out <- CP_Inference(
#'   h = h,
#'   i = 1,
#'   j = 2,
#'   Y = Y,
#'   res.CP.DPI = res.CP.DPI
#' )
#'
#' out$aij.h.de
#' out$se.h.ij
#' out$aij.iter
#' out$vartheta.ij
#' }
#'
#' @export
CP_Inference = function(h,i,j,Y,res.CP.DPI){
  
  A = res.CP.DPI$res.iter$A.hat
  Sigma.yij.xii.1 = res.CP.DPI$res.iter$Sigma.yij.xii.1
  f = res.CP.DPI$res.iter$f_hat
  n = NROW(f)
  
  aij.debias = aij.debias.iter(A,i,j,Sigma.yij.xii.1)
  
  aij = A[[j]][,i]
  
  aij.de = aij.debias$aij.de
  vartheta.ij = aij.debias$vartheta_ij
  
  se.ij  =  sqrt(cov.aij.debias.iter.est(h,i,j,A,f,Y)/n)
  
  aij.h.de = as.numeric(t(h)%*%aij.de)
  
  return(list(aij.h.de = aij.h.de,  se.h.ij = se.ij,  aij.iter = aij, vartheta.ij = vartheta.ij))
}



cp_residuals_general <- function(Y, f, A_list) {
  # Y: array of dimension n x d1 x d2 x ... x dm
  # f: matrix of dimension n x r
  # A_list: list(A1, ..., Am), where Aj is dj x r
  
  if (!is.array(Y)) {
    stop("Y must be an array with dimensions n x d1 x ... x dm.")
  }
  if (!is.matrix(f)) {
    stop("f must be an n x r matrix.")
  }
  if (!is.list(A_list) || length(A_list) < 1) {
    stop("A_list must be a non-empty list: list(A1, ..., Am).")
  }
  
  dims <- dim(Y)
  if (length(dims) < 2) {
    stop("Y must have at least 2 dimensions: n x d1.")
  }
  
  n <- dims[1]
  mode_dims <- dims[-1]
  m <- length(mode_dims)
  
  if (length(A_list) != m) {
    stop("length(A_list) must equal length(dim(Y)) - 1.")
  }
  
  r <- ncol(f)
  if (nrow(f) != n) {
    stop("nrow(f) must equal dim(Y)[1].")
  }
  
  # Check each loading matrix
  for (j in seq_len(m)) {
    Aj <- A_list[[j]]
    if (!is.matrix(Aj)) {
      stop(sprintf("A_list[[%d]] must be a matrix.", j))
    }
    if (nrow(Aj) != mode_dims[j]) {
      stop(sprintf("nrow(A_list[[%d]]) must equal dim(Y)[%d].", j, j + 1))
    }
    if (ncol(Aj) != r) {
      stop(sprintf("ncol(A_list[[%d]]) must equal ncol(f).", j))
    }
  }
  
  # Precompute rank-1 basis tensors (vectorized)
  total_dim <- prod(mode_dims)
  basis_mat <- matrix(0, nrow = total_dim, ncol = r)
  
  for (i in seq_len(r)) {
    # Start from the first mode loading vector
    comp_vec <- A_list[[1]][, i]
    
    # Sequentially build the tensor product
    if (m >= 2) {
      for (j in 2:m) {
        comp_vec <- as.vector(outer(comp_vec, A_list[[j]][, i]))
      }
    }
    
    basis_mat[, i] <- comp_vec
  }
  
  # Fitted values in matricized form: n x (d1*...*dm)
  Y_hat_mat <- f %*% t(basis_mat)
  
  # Convert back to array
  Y_hat <- array(Y_hat_mat, dim = dims)
  E <- Y - Y_hat
  
  list(
    residual = E,
    fitted = Y_hat,
    basis = basis_mat
  )
}
svd_inverse <- function(mat, threshold = 1e-6) {
  # Input validation
  if (!is.matrix(mat)) {
    stop("Input must be a matrix!")
  }
  
  if (!is.numeric(mat)) {
    stop("Matrix must be a numeric matrix!")
  }
  
  if (nrow(mat) != ncol(mat)) {
    stop("Only square matrices are supported for inversion!")
  }
  
  if (!is.numeric(threshold) || threshold <= 0) {
    stop("Threshold must be a positive number!")
  }
  
  # Perform singular value decomposition (SVD)
  svd_result <- svd(mat)
  
  # Extract singular value vector
  d <- svd_result$d
  
  # Replace singular values smaller than threshold with the threshold
  d_corrected <- pmax(d, threshold)
  
  # Construct the inverse of the diagonal matrix with corrected singular values
  d_inv <- diag(1 / d_corrected, nrow = length(d_corrected))
  
  # Calculate the corrected inverse matrix (V * D^{-1} * U^T)
  inv_mat <- svd_result$v %*% d_inv %*% t(svd_result$u)
  
  return(inv_mat)
}


# Function: ensure the first element of each column is positive
make_first_row_positive <- function(mat) {
  stopifnot(is.matrix(mat))  # ensure input is a matrix
  
  for (j in seq_len(ncol(mat))) {
    if (mat[1, j] < 0) {
      mat[, j] <- -mat[, j]
    }
  }
  return(mat)
}


generate_tensor_ar1_indep <- function(n, dims, c,
                                      error_dist = c("normal", "t"),
                                      df = NULL,
                                      burnin = 200,
                                      seed = NULL,
                                      standardize_t = TRUE,
                                      return_phi = TRUE) {
  # Basic checks
  if (!is.null(seed)) set.seed(seed)
  
  error_dist <- match.arg(error_dist)
  
  if (length(n) != 1 || !is.numeric(n) || n <= 0 || n != as.integer(n)) {
    stop("n must be one positive integer.")
  }
  
  if (length(dims) < 1 || any(!is.numeric(dims)) || any(dims <= 0) ||
      any(dims != as.integer(dims))) {
    stop("dims must be a vector of positive integers.")
  }
  
  if (length(c) != 1 || !is.numeric(c) || c < 0 || c >= 1) {
    stop("c must be a number in [0, 1).")
  }
  
  if (error_dist == "t") {
    if (is.null(df) || !is.numeric(df) || length(df) != 1 || df <= 0) {
      stop("For t errors, df must be one positive number.")
    }
    if (standardize_t && df <= 2) {
      stop("If standardize_t = TRUE, df must be greater than 2.")
    }
  }
  
  # Number of spatial locations
  p <- prod(dims)
  total_n <- n + burnin
  
  # Draw one AR coefficient for each location
  phi <- runif(p, min = -c, max = c)
  
  # Generate innovations matrix: total_n x p
  if (error_dist == "normal") {
    eps <- matrix(rnorm(total_n * p), nrow = total_n, ncol = p)
  } else {
    eps <- matrix(stats::rt(total_n * p, df = df), nrow = total_n, ncol = p)
    if (standardize_t) {
      # Scale t innovations to have variance 1
      eps <- eps / sqrt(df / (df - 2))
    }
  }
  
  # Allocate matrix for all AR(1) paths
  # Each column is one location-specific AR(1) process
  Y_mat <- matrix(0, nrow = total_n, ncol = p)
  
  # Initialize
  Y_mat[1, ] <- eps[1, ]
  
  # Time recursion, vectorized over all locations
  for (tt in 2:total_n) {
    Y_mat[tt, ] <- phi * Y_mat[tt - 1, ] + eps[tt, ]
  }
  
  # Remove burn-in
  Y_mat <- Y_mat[(burnin + 1):total_n, , drop = FALSE]
  
  # Reshape to tensor: c(n, d1, ..., dm)
  Y_tensor <- array(Y_mat, dim = c(n, dims))
  
  if (return_phi) {
    phi_tensor <- array(phi, dim = dims)
    return(list(
      Y = Y_tensor,      # Tensor time series
      phi = phi_tensor   # AR coefficients at each location
    ))
  } else {
    return(Y_tensor)
  }
}


Mat.k = function (A, k, eps = 10^-6)
{
  ev = eigen(A)$values
  mark = which(ev > eps)
  ev = ev[mark]
  evc = as.matrix(eigen(A)$vectors)[, mark]
  Matk = evc %*% diag(ev^k) %*% t(evc)
  return(Matk)
}






rho2.f.loss = function(f_hat,f){
  
  max( apply(1- cor(f_hat,f)^2 ,2,min))
  
}


rho2.loss.list = function(A_hat,A){
  
  m = length(A_hat)
  
  rho2  = c()
  for (j in 1:m) {
    rho2[j] = max(apply(1-(t(A_hat[[j]])%*%A[[j]])^2,2,min))
  }
  return(rho2)
}


vecpsi.loss = function(A_hat,A){
  apply(1-(t(A_hat)%*%A)^2,2,min)
}

vecpsi.loss.list = function(A_hat,A){
  m = length(A_hat)
  rho2  = vector()
  for (j in 1:m) {
    rho2  = rbind(rho2, vecpsi.loss(A_hat[[j]],A[[j]]))
  }
  apply(rho2, 2, max)
}


DGP.TCP <- function(n, m, D, r, w, ar.coef,
                    factor.loading = c("sparse-random", "sparse-random-corr1", "sparse-random-corr2"),
                    factor.corr = 0,
                    alpha = 0,
                    delta = 0.25,
                    par_E = 1,
                    heavytail = FALSE,
                    error.ar = FALSE) {
  factor.loading <- match.arg(factor.loading)
  
  if (alpha > 0) {
    if (factor.loading == "sparse-random") {
      A <- list()
      for (j in 1:m) {
        tau <- 0
        while (0 %in% tau || 1 %in% tau) {
          Aj <- matrix(runif(r * D[j], -1, 1), D[j], r)
          Aj[which(abs(Aj) < alpha * 1)] <- 0
          tau <- apply(Aj, 2, function(x) length(which(x != 0)))
        }
        A[[j]] <- apply(Aj, 2, l2s)
      }
    }
    
    if (factor.loading == "sparse-random-corr1") {
      A <- list()
      for (j in 1:m) {
        Aj <- matrix(runif(r * D[j], -1, 1), D[j], r)
        for (i in 2:r) {
          Aj[, i] <- delta * Aj[, i - 1] + Aj[, i]
        }
        for (i in 1:r) {
          Aj[sample(D[j], alpha * D[j]), i] <- 0
        }
        A[[j]] <- apply(Aj, 2, l2s)
      }
    }
    
    if (factor.loading == "sparse-random-corr2") {
      A <- list()
      for (j in 1:m) {
        Aj <- matrix(runif(r * D[j], -1, 1), D[j], r)
        Aj[sample(D[j], alpha * D[j]), ] <- 0
        for (i in 2:r) {
          Aj[, i] <- delta * Aj[, i - 1] + Aj[, i]
        }
        A[[j]] <- apply(Aj, 2, l2s)
      }
    }
  } else {
    A <- list()
    for (j in 1:m) {
      Aj <- matrix(runif(r * D[j], -3, 3), D[j], r)
      for (i in 2:r) {
        Aj[, i] <- delta * Aj[, i - 1] + Aj[, i]
      }
      A[[j]] <- apply(Aj, 2, l2s)
    }
  }
  
  f_m <- matrix(NA, n, r)
  S <- 0
  for (ii in 1:r) {
    xd <- arima.sim(model = list(ar = ar.coef[[ii]]), n = n)
    f_m[, ii] <- xd
  }
  
  COV <- diag(1 - factor.corr, r, r) + matrix(factor.corr, r, r)
  COV.half <- eigen(COV)$vectors %*% diag(sqrt(eigen(COV)$values)) %*% t(eigen(COV)$vectors)
  f_m <- f_m %*% COV.half
  
  for (ii in 1:r) {
    S_r <- f_m[, ii]
    for (jj in 1:m) {
      S_r <- S_r %o% A[[jj]][, ii]
    }
    S <- S + w[ii] * S_r
  }
  
  if (isFALSE(error.ar)) {
    if (isFALSE(heavytail)) {
      E <- array(rnorm(prod(D) * n, 0, par_E), c(n, D))
    } else {
      E <- array(stats::rt(prod(D) * n, heavytail), c(n, D))
    }
  } else {
    if (isFALSE(heavytail)) {
      resE <- generate_tensor_ar1_indep(n = n, dims = D, c = error.ar, error_dist = "normal")
    } else {
      resE <- generate_tensor_ar1_indep(n = n, dims = D, c = error.ar, error_dist = "t", df = 5)
    }
    E <- resE$Y
  }
  
  list(Y = Y <- S + E, C = S, E = E, A = A, f = f_m, w = w)
}


Mat.tensor = function(Y,j){ # fold on j-th mode
  t(apply(Y,j,c))
}

Threshold.Tensor = function(SigmaY,n,sigma0,delta){
  d = prod(dim(SigmaY))
  SigmaY[which(abs(SigmaY) < delta*sigma0*sqrt(log(d)/n))] <- 0
  return(SigmaY)
}



Autocov_xi_Y_nothres = function(Y,eta,lag.k = k){
  
  if(length(dim(Y)) == 3){
    n = dim(Y)[1]
    k = lag.k
    
    Y_mean = 0
    eta_mean = 0
    for (ii in 1:n) {
      Y_mean = Y_mean + Y[ii,,]
      eta_mean = eta_mean + eta[ii]
    }
    Y_mean   = Y_mean/n
    eta_mean = eta_mean/n
    
    Sigma_Y_eta_k = 0
    for (ii in (k+1):n) {
      Sigma_Y_eta_k = Sigma_Y_eta_k + (Y[ii,,] - Y_mean)*(eta[ii-k] - eta_mean)
    }
    
  }
  
  if(length(dim(Y)) == 2){
    n = dim(Y)[1]
    k = lag.k
    
    Y_mean = 0
    eta_mean = 0
    for (ii in 1:n) {
      Y_mean = Y_mean + Y[ii,]
      eta_mean = eta_mean + eta[ii]
    }
    Y_mean   = Y_mean/n
    eta_mean = eta_mean/n
    
    Sigma_Y_eta_k = 0
    for (ii in (k+1):n) {
      Sigma_Y_eta_k = Sigma_Y_eta_k + (Y[ii,] - Y_mean)*(eta[ii-k] - eta_mean)
    }
    
  }
  
  return(Sigma_Y_eta_k/(n-k))
}



tensor.Autocov_xi_Y = function(Y,xi,lag.k = k){ # reuturn a autocovariance  tensor (d1 x d2 x ... x d_m)
  
  n = dim(Y)[1]
  k = lag.k
  
  Y_mean = apply(Y, c(2:length(dim(Y))), mean)
  xi_mean = mean(xi)
  
  Sigma_Y_xi_k = 0
  for (ii in (k+1):n) {
    Sigma_Y_xi_k = Sigma_Y_xi_k + (base_extract(Y, 1, ii, drop = TRUE) - Y_mean)*(xi[ii-k] - xi_mean)
  }
  return(Sigma_Y_xi_k/(n-k))
}


HDTTS.CP.est =  function(Y,
                         xi = NULL,
                         Rank = NULL,
                         K = 10,
                         Threshold = FALSE,
                         delta = NULL,
                         Ratio.type ="log",
                         grid_delta1 = 50,
                         delta_max = 0.1){
  
  est.ABr =   function(Sigma.tensor.Y.k_list_orginal,Rank,delta,sigma0){
    Sigma.tensor.Y.k_list = list()
    
    for (kk in 1:K) {
      Sigma.tensor.Y.k = Sigma.tensor.Y.k_list_orginal[[kk]]
      
      if(Threshold == TRUE & kk == 1){
        Sigma.tensor.Y.k = Threshold.Tensor(Sigma.tensor.Y.k,n = n,sigma0 = sigma0,delta = delta)
        index = Sigma.tensor.Y.k
      }
      if(Threshold == TRUE & kk > 1){
        Sigma.tensor.Y.k[which(index == 0)] <- 0
      }
      
      Sigma.tensor.Y.k_list[[kk]] = Sigma.tensor.Y.k
    }
    
    r_tol  = vector()
    Mjk_list = MMjk_list = list()
    Sigma.Y.k_list_tol = list()
    for (j in 1:m){
      Mjk = 0
      MMjk = 0
      Sigma.Y.k_list = list()
      for (kk in 1:K){
        Sigma.Y.k = Mat.tensor(Sigma.tensor.Y.k_list[[kk]],j)
        Mjk = Mjk + Sigma.Y.k%*%t(Sigma.Y.k)
        MMjk = MMjk + t(Sigma.Y.k)%*%Sigma.Y.k
        Sigma.Y.k_list[[kk]] = Sigma.Y.k
      }
      
      Mjk_list[[j]] = Mjk
      MMjk_list[[j]] = MMjk
      Sigma.Y.k_list_tol[[j]] = Sigma.Y.k_list
      
      if(Ratio.type == "log"){
        eigenvalue_j =  log(eigen(Mjk)$values + 1)
      }
      if(Ratio.type == "classical"){
        eigenvalue_j =  eigen(Mjk)$values
      }
      
      ratio_j = (eigenvalue_j[-1] + sigma0/n)/(eigenvalue_j[-length(eigenvalue_j)] + sigma0/n)
      
      r_j = which.min(ratio_j[1:(0.5*length(ratio_j))])
      
      r_tol = c(r_tol,r_j)
    }
    
    length_enc <- rle(r_tol)
    
    r =  max(r_tol)
    r.hat = max(r_tol)
    
    if(!is.null(Rank)){
      r = Rank
    }
    
    A_hat = A12_hat = A1K_hat = list()
    j = 1
    
    P_list = list()
    Q_list = list()
    K_12_tilde_list = list()
    
    for (j in 1:m) {
      Pj = as.matrix(eigen(Mjk_list[[j]])$vectors[,1:r])
      Qj = as.matrix(eigen(MMjk_list[[j]])$vectors[,1:r])
      
      P_list[[j]] = Pj
      Q_list[[j]] = Qj
      
      Sigma.Y.k_list = Sigma.Y.k_list_tol[[j]]
      
      
      bb1 = Sigma.Y.k_list[[1]]%*%Qj
      bb2 = Sigma.Y.k_list[[2]]%*%Qj
      bbk = Sigma.Y.k_list[[K]]%*%Qj
      
      K21_j_tilde = bb2%*%MASS::ginv(t(bb1)%*%bb1)%*%t(bb1)
      
      K12_j_tilde = bb1%*%MASS::ginv(t(bb2)%*%bb2)%*%t(bb2)
      
      K1K_j_tilde = bb1%*%MASS::ginv(t(bbk)%*%bbk)%*%t(bbk)
      
      evd.K12 = eigen(K12_j_tilde)
      
      if(sum(abs(Im(evd.K12$values)[1:r])) < 10^-10){
        Aj_hat = Re(eigen(K21_j_tilde)$vectors[,1:r])
        
        Aj12_hat = Re(eigen(K12_j_tilde)$vectors[,1:r])
        
        Aj1K_hat = Re(eigen(K1K_j_tilde)$vectors[,1:r])
      }else{
        Aj_hat = Aj1K_hat = Aj12_hat =  apply(Complex2Real(eigen(K12_j_tilde)$vectors[,1:r]) ,2, l2s)
      }
      
      
      
      K_12_tilde_list[[j]] = K12_j_tilde
      
      A_hat[[j]] = as.matrix(Aj_hat)
      
      A12_hat[[j]] = as.matrix(Aj12_hat)
      
      A1K_hat[[j]] = as.matrix(Aj1K_hat)
    }
    
    delta_sel = delta
    
    return(list(A.hat = A12_hat, r.hat = r.hat, delta_sel = delta,Sigma.tensor.Y.k_list_threshold = Sigma.tensor.Y.k_list,Q_list = Q_list, K_12_tilde_list = K_12_tilde_list))
  }
  
  
  n = dim(Y)[1]
  D = dim(Y)[-1]
  m = length(D)
  
  if(is.null(xi)){
    xi = tensor.est.xi(Y,random = F)
  }
  
  sigma0 = sqrt(sum(Y^2)/(n*prod(D)))
  
  Sigma.tensor.Y.k_list_orginal = list()
  
  for (kk in 1:K) {
    Sigma.tensor.Y.k_list_orginal[[kk]]  = tensor.Autocov_xi_Y(Y,scale(xi),kk)
  }
  
  if(is.null(delta) & Threshold == TRUE){
    
    test_list = z_list = r_list  = vector()
    
    net_delta = seq(0,delta_max,length.out = grid_delta1)
    
    A_hat_list = list()
    
    for (delta in net_delta) {
      
      Sigma.tensor.Y.k_list = list()
      
      for (kk in 1:K) {
        Sigma.tensor.Y.k = Sigma.tensor.Y.k_list_orginal[[kk]]
        
        if(Threshold == TRUE & kk == 1){
          Sigma.tensor.Y.k = Threshold.Tensor(Sigma.tensor.Y.k,n = n,sigma0 = sigma0,delta = delta)
          index = Sigma.tensor.Y.k
        }
        if(Threshold == TRUE & kk > 1){
          Sigma.tensor.Y.k[which(index == 0)] <- 0
        }
        Sigma.tensor.Y.k_list[[kk]] = Sigma.tensor.Y.k
      }
      
      test_tol = z_tol = r_tol  = vector()
      
      A_hat = list()
      
      for (j in 1:m) {
        Mjk = 0
        MMjk = 0
        Sigma.Y.k_list = list()
        for (kk in 1:K){
          Sigma.Y.k = Mat.tensor(Sigma.tensor.Y.k_list[[kk]],j)
          Mjk = Mjk + Sigma.Y.k%*%t(Sigma.Y.k)
          MMjk = MMjk + t(Sigma.Y.k)%*%Sigma.Y.k
          Sigma.Y.k_list[[kk]] = Sigma.Y.k
        }
        
        if(Ratio.type == "log"){
          eigenvalue_j =  log(eigen(Mjk)$values + 1)
        }
        if(Ratio.type == "classical"){
          eigenvalue_j =  eigen(Mjk)$values
        }
        
        ratio_j = (eigenvalue_j[-1] + sigma0/n)/(eigenvalue_j[-length(eigenvalue_j)] + sigma0/n)
        
        ratio_j = ratio_j[1:(0.5*length(ratio_j))]
        
        r_j = which.min(ratio_j)
        
        test_j =   eigenvalue_j[r_j]
        
        z_j = min(ratio_j)
        
        test_tol = c(test_tol,test_j)
        z_tol = c(z_tol,z_j)
        
        r_tol = c(r_tol,r_j)
      }
      
      test_list = rbind(test_list, test_tol)
      z_list = rbind(z_list, z_tol)
      r_list = rbind(r_list, r_tol)
      
    }
    
    place1 = which.min(rowMeans(z_list))
    
    delta_sel_1 = net_delta[place1]
    
    
    res.1 = est.ABr(Sigma.tensor.Y.k_list_orginal,Rank = Rank,delta = delta_sel_1,sigma0)
    
    
  }else{ #no thresholding or given delta thresholding
    z_list = r_list =  test_list = NULL
    
    res.1 = est.ABr(Sigma.tensor.Y.k_list_orginal,Rank = Rank,delta = delta,sigma0)
  }
  
  res.1$Sigma.tensor.Y.k_list_orginal = Sigma.tensor.Y.k_list_orginal
  res.1$sigma0 = sigma0
  
  return( res.1 )
  
}



aij.debias.iter  = function(A, i,j, Sigma.yij.xii.1){
  m = length(A)
  r = NCOL(A[[1]])
  aij = A[[j]][,i]
  
  vartheta_ij = (as.numeric(t(aij)%*% Sigma.yij.xii.1[[j]][,i])*aij - Sigma.yij.xii.1[[j]][,i])/as.numeric(t(aij)%*% Sigma.yij.xii.1[[j]][,i])
  
  
  aij.de = aij - vartheta_ij
  
  return(list(aij.de = aij.de, vartheta_ij = vartheta_ij))
}




cov.aij.debias.iter.est = function(h,i,j,A,f,Y){
  aij = A[[j]][,i]
  
  m  = length(A)
  n  = NROW(f)
  dj = length(aij)
  
  Aj = A[[j]]
  
  if(m >= 3){
    Bj = rTensor::khatri_rao_list(A[c(m:1)[-(m-j+1)]])
  }else{
    Bj = A[c(1:m)[-j]][[1]]
  }
  
  B.MP.j = Bj%*%MASS::ginv(t(Bj)%*%Bj)
  
  ff = scale(f)
  
  xi.vmax.j = lm(ff[-n,i] ~ ff[-1,-i]-1)$residuals
  
  sigma.fi.xi.wi = Autocov_xi_Y_nothres(f,c(xi.vmax.j,0),1)[i]
  
  beta = (B.MP.j[,i]) %x% t(t(h)%*%(diag(dj) - aij%*%t(aij)))
  
  Wj = rTensor::khatri_rao(Bj,Aj)
  
  Wj.MP = Wj%*%MASS::ginv(t(Wj)%*%Wj)
  
  Yp <- aperm(Y, c(1,j+1,c(2:(m+1))[-j])) # change the mode order into n x dj x d1 x ... x dm
  
  Yj = Mat.tensor(Y,1)
  
  Ej = Yj - Yj%*%Wj.MP%*%t(Wj)
  
  qj  = Ej%*%beta
  qj2 = Yj%*%beta
  
  cov.tol = mean(xi.vmax.j^2*(qj[2:n])^2)
  
  COV.TOL = cov.tol/sigma.fi.xi.wi^2
  
  return(as.numeric(COV.TOL))
}


CP.iter.DPI.xi <- function(A.hat,
                           K,
                           Y,
                           n,
                           delta2,
                           sigma0,
                           iter_max = 20,
                           eps = 1e-5,
                           print.eps = TRUE,
                           A = NULL,
                           iter_lag = 1,
                           all.put = FALSE,
                           Component = NULL
) {
  
  m <- length(A.hat)
  D <- dim(Y)[-1]
  r <- NCOL(A.hat[[1]])
  
  A.tol <- A.1 <- A.0 <- Sigma.yij.xii.1 <- A.hat
  A.tol.list <- list()
  wps.tol <- vector()
  iter.error.mat <- matrix(0, iter_max + 1, 1)
  fnorm.resid <- vector()
  Uj.mat.inl <- NULL
  A.tilde.j.inl <- NULL
  B.tilde.j.inl <- NULL
  Yp_hat <- NULL
  
  for (ll in 1:iter_max) {
    if (!is.null(A)) {
      iter.error.mat[c(ll:(iter_max + 1)), ] <- matrix(
        max(rho2.loss.list(A.0, A)),
        length(c(ll:(iter_max + 1))),
        1,
        byrow = TRUE
      )
    }
    
    A.tol <- A.0
    
    for (j in 1:m) {
      delta_2j <- delta2[j]
      Yp <- aperm(Y, c(1, j + 1, c(2:(m + 1))[-j]))
      
      if (m >= 3) {
        B.tilde.j <- rTensor::khatri_rao_list(A.0[c(m:1)[-(m - j + 1)]])
      } else {
        B.tilde.j <- A.0[c(1:m)[-j]][[1]]
      }
      A.tilde.j <- A.0[[j]]
      
      B.MP.j <- B.tilde.j %*% MASS::ginv(t(B.tilde.j) %*% B.tilde.j)
      A.MP.j <- A.tilde.j %*% MASS::ginv(t(A.tilde.j) %*% A.tilde.j)
      
      if (m >= 3) {
        dim(Yp) <- c(dim(Yp)[1:2], prod(dim(Yp)[(m + 1):3]))
      }
      
      Uj <- rTensor::ttl(rTensor::as.tensor(Yp), list_mat = list(t(A.MP.j), t(B.MP.j)), ms = c(2, 3))@data
      
      if (NCOL(Uj) > 1) {
        Uj.mat <- t(apply(Uj, 1, diag))
      } else {
        Uj.mat <- as.matrix(Uj)
      }
      
      if (j == 1) {
        Yp_hat.tmp <- 0
        for (k in 1:r) {
          Yp_hat.tmp <- Yp_hat.tmp + Uj.mat[, k] %o% A.tilde.j[, k] %o% B.tilde.j[, k]
        }
        fnorm.resid[ll] <- fnorm(Yp - Yp_hat.tmp) / fnorm(Yp)
        
        if (ll == 1) {
          Uj.mat.inl <- as.matrix(Uj.mat)
          A.tilde.j.inl <- A.tilde.j
          B.tilde.j.inl <- B.tilde.j
        }
      }
      
      Uj.reg <- as.matrix(Uj.mat)
      
      for (k in 1:r) {
        A.fnorm <- vector()
        A.breve.jk.list <- vector()
        
        for (qq in 1:iter_lag) {
          if (r > 1) {
            ujk <- lm(Uj.reg[-c((n - qq + 1):n), k] ~ Uj.reg[-c(1:qq), -k] - 1)$residuals
            sigma.k.ujk <- tensor.Autocov_xi_Y(Y, c(scale(ujk), rep(0, qq)), qq)
          } else {
            ujk <- Uj.reg
            sigma.k.ujk <- tensor.Autocov_xi_Y(Y, scale(ujk), qq)
          }
          
          A.breve.jk <- Mat.tensor(sigma.k.ujk, j) %*% B.MP.j[, k]
          A.fnorm[qq] <- fnorm(A.breve.jk)
          A.breve.jk.list <- as.matrix(cbind(A.breve.jk.list, A.breve.jk))
        }
        
        A.breve.jk <- as.matrix(A.breve.jk.list[, which.max(A.fnorm)])
        
        Sigma.yij.xii.1[[j]][, k] <- A.breve.jk
        
        A.breve.jk.thres <- Threshold.Tensor(A.breve.jk, n, sigma0, delta_2j)
        
        if (sum(A.breve.jk.thres) < 0.01) {
          A.breve.jk.thres <- Threshold.Tensor(A.breve.jk, n, sigma0, 0)
        }
        
        
        A.1[[j]][, k] <- apply(A.breve.jk.thres, 2, l2s)
      }
      
      A.0 <- A.1
    }
    
    A.tol.new <- A.1
    A.tol.list[[ll]] <- A.tol.new
    
    wps <- sum(sqrt(abs(rho2.loss.list(A.tol.new, A.tol))))
    wps.tol <- rbind(wps.tol, wps)
    
    if (isTRUE(print.eps)) {
      cat("\n", round(wps, 6))
    }
    
    if (wps < eps) {
      break
    }
  }
  
  iter.error <- c(iter.error.mat)
  A.tol <- A.tol.new
  
  if (ll == iter_max) {
    A.tol <- A.tol.list[[which.min(fnorm.resid)]]
  }
  
  j <- 1
  Yp <- aperm(Y, c(1, j + 1, c(2:(m + 1))[-j]))
  if (m >= 3) {
    B.tilde.j <- rTensor::khatri_rao_list(A.tol[c(m:1)[-(m - j + 1)]])
  } else {
    B.tilde.j <- A.tol[c(1:m)[-j]][[1]]
  }
  A.tilde.j <- A.tol[[j]]
  
  B.MP.j <- B.tilde.j %*% MASS::ginv(t(B.tilde.j) %*% B.tilde.j)
  A.MP.j <- A.tilde.j %*% MASS::ginv(t(A.tilde.j) %*% A.tilde.j)
  if (m >= 3) {
    dim(Yp) <- c(dim(Yp)[1:2], prod(dim(Yp)[(m + 1):3]))
  }
  
  Uj <- rTensor::ttl(rTensor::as.tensor(Yp), list_mat = list(t(A.MP.j), t(B.MP.j)), ms = c(2, 3))@data
  
  if (r > 1) {
    Uj.mat <- t(apply(Uj, 1, diag))
  } else {
    Uj.mat <- as.matrix(Uj)
  }
  
  CP_loss_iter <- CP_loss_inl <- -999
  if (!is.null(Component)) {
    Cp <- aperm(Component, c(1, j + 1, c(2:(m + 1))[-j]))
    Yp_hat <- 0
    Yp_hat_inl <- 0
    for (k in 1:r) {
      Yp_hat <- Yp_hat + Uj.mat[, k] %o% A.tilde.j[, k] %o% B.tilde.j[, k]
      Yp_hat_inl <- Yp_hat_inl + Uj.mat.inl[, k] %o% A.tilde.j.inl[, k] %o% B.tilde.j.inl[, k]
    }
    CP_loss_iter <- fnorm(Yp_hat - Cp) / sqrt(n * prod(D))
    CP_loss_inl <- fnorm(Yp_hat_inl - Cp) / sqrt(n * prod(D))
  }
  
  if (isFALSE(all.put)) {
    res.list <- list(
      A.hat = A.tol,
      A.init = A.hat,
      Sigma.yij.xii.1 = Sigma.yij.xii.1,
      r = r,
      iter_step = ll,
      fnorm.resid = fnorm.resid,
      f_hat = as.matrix(Uj.mat),
      f_hat_inl = as.matrix(Uj.mat.inl)
    )
  } else {
    res.list <- list(
      A.hat = A.tol,
      A.init = A.hat,
      A.tol.list = A.tol.list,
      Sigma.yij.xii.1 = Sigma.yij.xii.1,
      r = r,
      delta2_sel = delta2,
      iter_step = ll,
      iter_error = iter.error,
      wps.tol = wps.tol,
      fnorm.resid = fnorm.resid,
      f_hat = Uj.mat,
      f_hat_inl = Uj.mat.inl,
      CP_loss = c(CP_loss_inl, CP_loss_iter),
      Yp_hat = Yp_hat
    )
  }
  
  res.list
}


cov.aij.debias.iter <- function(h, A, i, j, f, COV.VECE, w) {
  m <- length(A)
  D <- sapply(A, nrow)
  aij <- A[[j]][, i]
  dj <- length(aij)
  Aj <- A[[j]]
  
  if (m >= 3) {
    B.tilde.j <- rTensor::khatri_rao_list(A[c(m:1)[-(m - j + 1)]])
  } else {
    B.tilde.j <- A[c(1:m)[-j]][[1]]
  }
  
  B.MP.j <- B.tilde.j %*% MASS::ginv(t(B.tilde.j) %*% B.tilde.j)
  
  Fjk <- lm(f[-NROW(f), i] ~ f[-1, -i] - 1)$residuals
  sigma.fi.xi <- Autocov_xi_Y_nothres(f, c(Fjk, 0), 1)[i] / sqrt(var(Fjk))
  
  G <- t(h) %*% (diag(dj) - aij %*% t(aij)) %*% (t(B.MP.j[, i]) %x% diag(D[j]))
  COV.TOL <- G %*% COV.VECE %*% t(G) / sigma.fi.xi^2 / (w[i])^2
  
  as.numeric(COV.TOL)
}



tensor.est.xi = function(Y,d_max = 10,thresh_per = 0.99, random = F, seed = NULL){
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  n = dim(Y)[1];D = dim(Y)[-1]
  
  Y.mat = Mat.tensor(Y,1)
  
  if(n > prod(D)){
    eig_Y.mat = eigen(MatMult(t(Y.mat),Y.mat))
    cfr =  cumsum(eig_Y.mat$values)/sum(eig_Y.mat$values)
    d_hat = min(which(cfr > thresh_per))
    d_fin = min(d_max,d_hat)
    
    w_inl = as.matrix(eig_Y.mat$vectors[,1:d_fin])
    
    if(random == T){
      w_hat = (w_inl)%*%randortho(d_fin)
    }else{
      w_hat = (w_inl)
    }
    
    w_hat = make_first_row_positive(w_hat)
    
    xi.f = scale(Y.mat%*%w_hat)
    xi   = rowMeans(xi.f)
    
  }else{
    eig_Y.mat = eigen(MatMult(Y.mat,t(Y.mat)))
    cfr = cumsum(eig_Y.mat$values)/sum(eig_Y.mat$values)
    d_hat = min(which(cfr > thresh_per))
    d_fin = min(d_max,d_hat)
    
    w_inl = as.matrix(eig_Y.mat$vectors[,1:d_fin])
    
    if(random == T){
      xi.f = (w_inl)%*%randortho(d_fin)
    }else{
      xi.f = as.matrix(w_inl)
    }
    
    weight = sqrt(eig_Y.mat$values[1:d_fin])
    
    xi.f = xi.f   #%*%diag(weight)
    
    xi.f = make_first_row_positive(xi.f)
    
    xi = rowMeans(xi.f)
    
  }
  return(scale(xi))
}

RP.xi.sel <- function(Y,
                      r_breve = NULL,
                      eps = 0.1,
                      lag.k = 10,
                      Randomized.time = 50,
                      A = NULL) {
  n <- dim(Y)[1]
  D <- dim(Y)[-1]
  m <- length(D)
  
  Rank <- r_breve
  
  gg <- 0
  iter <- 1
  while (gg == 0 && iter < 20) {
    A.hat.NT_list <- list()
    xi_list <- list()
    
    for (ss in 1:Randomized.time) {
      xi <- tensor.est.xi(Y, random = TRUE)
      res.init <- HDTTS.CP.est(
        Y = Y,
        xi = xi,
        Rank = Rank,
        K = lag.k,
        Threshold = FALSE,
        delta = 0,
        Ratio.type = "log"
      )
      xi_list[[ss]] <- xi
      A.hat.NT_list[[ss]] <- res.init$A.hat
    }
    
    G <- matrix(1, Randomized.time, Randomized.time)
    H <- vector()
    diag(G) <- 100
    eps0 <- eps
    DD_tol <- vector()
    
    for (vv in 1:r_breve) {
      G_i <- matrix(0, Randomized.time, Randomized.time)
      
      for (jj in 1:Randomized.time) {
        for (kk in 1:Randomized.time) {
          if (jj != kk) {
            loss_vec <- vecpsi.loss.list(A.hat.NT_list[[jj]], A.hat.NT_list[[kk]])
            H <- rbind(H, c(jj, kk, loss_vec))
            D_ijk <- abs(loss_vec[vv])
            G_i[jj, kk] <- ifelse(D_ijk < eps0, 1, 0)
          }
        }
      }
      
      DD_tol <- rbind(DD_tol, colSums(G_i))
    }
    
    gg <- max(colSums(DD_tol))
    iter <- iter + 1
  }
  
  if (iter == 20) {
    xi.sel <- tensor.est.xi(Y, d_max = 1)
  } else {
    place <- which.max(colSums(DD_tol))
    xi.sel <- xi_list[[place]]
  }
  
  list(xi.sel = xi.sel, H = H, G = G)
}

randortho <- function(n, type = c("orthonormal", "unitary")) {
  stopifnot(
    is.numeric(n),
    length(n) == 1,
    floor(n) == ceiling(n),
    n >= 1
  )
  
  if (n == 1) {
    return(matrix(1, 1, 1))
  }
  
  type <- match.arg(type)
  
  # Generate a random real or complex Gaussian matrix
  if (type == "orthonormal") {
    z <- matrix(rnorm(n * n), nrow = n, ncol = n) / sqrt(2)
  } else {
    z <- (matrix(rnorm(n * n), nrow = n, ncol = n) +
            1i * matrix(rnorm(n * n), nrow = n, ncol = n)) / sqrt(2)
  }
  
  # QR decomposition
  Z <- qr(z)
  q <- qr.Q(Z)
  r <- qr.R(Z)
  
  # Adjust the phases/signs to ensure uniformity
  d <- diag(r)
  ph <- d / abs(d)
  
  q %*% diag(ph)
}


base_extract <- function(Y, dim = 1, index, drop = TRUE) {
  args <- rep(list(substitute()), length(dim(Y)))
  args[[dim]] <- index
  args$drop <- drop
  do.call("[", c(list(Y), args))
}


