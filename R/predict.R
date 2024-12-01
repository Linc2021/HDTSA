#' Make predictions from a \code{"factors"} object
#'
#' This function makes predictions from a \code{"factors"} object.
#' 
#' @details
#' Forecasting for \eqn{{\bf y}_t} can be implemented in two steps:
#' 
#' \emph{Step 1}. Get the \eqn{h}-step ahead forecast of the \eqn{\hat{r} \times 1} 
#' time series \eqn{\hat{\bf x}_t} [See \code{\link{Factors}}], denoted by
#' \eqn{\hat{\bf x}_{n+h}}, using a VAR model 
#' (if \eqn{\hat{r} > 1}) or an ARIMA model (if \eqn{\hat{r} = 1}). The orders 
#' of VAR and ARIMA models are determined by AIC by default. Otherwise, they 
#'  can also be specified by users through the arguments \code{control_VAR}
#'  and \code{control_ARIMA}, respectively.
#' 
#' \emph{Step 2}. The forecasted value for \eqn{{\bf y}_t} is obtained by
#'  \eqn{\hat{\bf y}_{n+h}= \hat{\bf A}\hat{\bf x}_{n+h}}.
#' 
#' 
#'
#' @param object An object of class \code{"factors"} constructed by \code{\link{Factors}}.
#' @param newdata Optional. A new data matrix to predict from.
#' @param ... Currently not used.
#' @param control_ARIMA A list of arguments passed to the function
#' \code{auto.arima()} of \pkg{forecast}. See 'Details' and the manual of \code{auto.arima()}.
#' The default is \code{list(ic = "aic")}.
#' @param control_VAR A list of arguments passed to the function
#' \code{VAR()} of \pkg{vars}. See 'Details' and the manual of \code{VAR()}.
#' The default is \code{list(type = "const", lag.max = 6, ic = "AIC")}.
#' @param n.ahead An integer specifying the number of steps ahead for prediction.
#' 
#' @return \item{ts_pred}{A matrix of predicted values.}
#' @seealso \code{\link{Factors}}
#' @examples
#' library(HDTSA)
#' data(FamaFrench, package = "HDTSA")
#' 
#' ## Remove the market effects
#' reg <- lm(as.matrix(FamaFrench[, -c(1:2)]) ~ as.matrix(FamaFrench$MKT.RF))
#' Y_2d = reg$residuals
#' 
#' res_factors <- Factors(Y_2d, lag.k = 5)
#' pred_fac_Y <- predict(res_factors, n.ahead = 1)
#' 
#' @importFrom vars VAR
#' @importFrom forecast auto.arima
#' @importFrom stats predict
#' @method predict factors
#' @export
predict.factors <- function(object, newdata = NULL, n.ahead = 10,
                            control_ARIMA = list(), control_VAR = list(), ...) {
  params <- list(...)
  con1 <- list(ic = "aic")
  con1[(namc <- names(control_ARIMA))] <- control_ARIMA
  
  con2 <- list(type = "const", lag.max = 6, ic = "AIC")
  con2[(namc <- names(control_VAR))] <- control_VAR
  
  if (!is.null(object$loading.mat)) {
    loading <- object$loading.mat
  }
  if (missing(newdata) || is.null(newdata)) {
    ts <- object$X
  }
  else {ts <- newdata %*% as.matrix(loading)}
  
  if (!is.null(object$reg.coff.mat)){
    stop("Not yet implemented for HDSreg.")
  }
  ts_pred <- Forecast(ts, con1, con2, n.ahead)$f.forecast %*% t(loading)
  ts_pred
}

#' Make predictions from a \code{"tspca"} object
#'
#' This function makes predictions from a \code{"tspca"} object.
#' 
#' @details
#' Having obtained \eqn{\hat{\bf x}_t} using the \code{\link{PCA_TS}} function, which is
#' segmented into \eqn{q} uncorrelated subseries
#' \eqn{\hat{\bf x}_t^{(1)},\ldots,\hat{\bf x}_t^{(q)}}, the forecasting for \eqn{{\bf y}_t}
#' can be performed in two steps:
#' 
#' \emph{Step 1}. Get the \eqn{h}-step ahead forecast \eqn{\hat{\bf x}_{n+h}^{(j)}} \eqn{(j=1,\ldots,q)}
#'  by using a VAR model (if the dimension of \eqn{\hat{\bf x}_t^{(j)}} is larger than 1) 
#'  or an ARIMA model (if the dimension of \eqn{\hat{\bf x}_t^{(j)}} is 1). The orders 
#'  of VAR and ARIMA models are determined by AIC by default. Otherwise, they 
#'  can also be specified by users through the arguments \code{control_VAR} and \code{control_ARIMA}, respectively.
#' 
#' \emph{Step 2}. Let \eqn{\hat{\bf x}_{n+h} = (\{\hat{\bf x}_{n+h}^{(1)}\}',\ldots,\{\hat{\bf x}_{n+h}^{(q)}\}')'}.
#' The forecasted value for \eqn{{\bf y}_t} is obtained by
#'  \eqn{\hat{\bf y}_{n+h}= \hat{\bf B}^{-1}\hat{\bf x}_{n+h}}.
#'
#' @param object An object of class \code{"tspca"} constructed by \code{\link{PCA_TS}}.
#' @param newdata Optional. A new data matrix to predict from.
#' @param ... Currently not used.
#' @param control_ARIMA A list of arguments passed to the function
#' \code{auto.arima()} of \pkg{forecast}. See 'Details' and the manual of \code{auto.arima()}.
#' The default is \code{list(max.d = 0, max.q = 0, ic = "aic")}.
#' @param control_VAR A list of arguments passed to the function
#' \code{VAR()} of \pkg{vars}. See 'Details' and the manual of \code{VAR()}.
#' The default is \code{list(type = "const", lag.max = 6, ic = "AIC")}.
#' @param n.ahead An integer specifying the number of steps ahead for prediction.
#' 
#' @seealso \code{\link{PCA_TS}}
#' @examples
#' library(HDTSA)
#' data(FamaFrench, package = "HDTSA")
#' 
#' ## Remove the market effects
#' reg <- lm(as.matrix(FamaFrench[, -c(1:2)]) ~ as.matrix(FamaFrench$MKT.RF))
#' Y_2d = reg$residuals
#' 
#' res_pca <- PCA_TS(Y_2d, lag.k = 5, thresh = TRUE)
#' pred_pca_Y <- predict(res_pca, n.ahead = 1)
#' 
#' 
#' @return \item{Y_pred}{A matrix of predicted values.}
#' @method predict tspca
#' @export
predict.tspca <- function(object, newdata = NULL, n.ahead = 10,
                          control_ARIMA = list(), control_VAR = list(), ...) {
  params <- list(...)
  con1 <- list(max.d = 0, max.q = 0, ic = "aic")
  con1[(namc <- names(control_ARIMA))] <- control_ARIMA
  
  con2 <- list(type = "const", lag.max = 6, ic = "AIC")
  con2[(namc <- names(control_VAR))] <- control_VAR
  
  X <- object$X
  B <- object$B
  group <- object$Groups

  if (missing(newdata) || is.null(newdata)) {
    newdata <- X
  }
  p <- ncol(newdata)
  ts_pred <- matrix(0, n.ahead, p)
  for (ii in seq_len(ncol(group))) {
    ts_i <- newdata[ , group[, ii], drop = F]
    ts_i_pred <- Forecast(ts_i, con1, con2, n.ahead)$f.forecast
    ts_pred[1:n.ahead, group[, ii]] <- ts_i_pred
  }
  Y_pred <- ts_pred %*% t(solve(B))
  Y_pred
}

#' Make predictions from a \code{"mtscp"} object
#'
#' This function makes predictions from a \code{"mtscp"} object.
#' 
#' @details
#' Forecasting for \eqn{{\bf y}_t} can be implemented in two steps:
#' 
#' \emph{Step 1}. Get the \eqn{h}-step ahead forecast of the \eqn{\hat{d} \times 1} 
#' time series \eqn{\hat{\bf x}_t=(\hat{x}_{t,1},\ldots,\hat{x}_{t,\hat{d}})'}
#' [See \code{\link{CP_MTS}}], denoted by \eqn{\hat{\bf x}_{n+h}}, using a VAR model 
#' (if \eqn{\hat{d} > 1}) or an ARIMA model
#' (if \eqn{\hat{d} = 1}). The orders of VAR and ARIMA models are determined by
#' AIC by default. Otherwise, they can also be specified by users through the
#' arguments \code{control_VAR} and \code{control_ARIMA}, respectively.
#' 
#' \emph{Step 2}. The forecasted value for \eqn{{\bf Y}_t} is obtained by
#'  \eqn{\hat{\bf Y}_{n+h}= \hat{\bf A} \hat{\bf X}_{n+h} \hat{\bf B}'} with
#'  \eqn{ \hat{\bf X}_{n+h} = {\rm diag}(\hat{\bf x}_{n+h})}.
#'
#' @param object An object of class \code{"mtscp"} constructed by \code{\link{CP_MTS}}.
#' @param newdata Optional. A new data matrix to predict from.
#' @param ... Currently not used.
#' @param control_ARIMA A list of arguments passed to the function
#' \code{auto.arima()} of \pkg{forecast}. See 'Details' and the manual of \code{auto.arima()}.
#' The default is \code{list(ic = "aic")}.
#' @param control_VAR A list of arguments passed to the function
#' \code{VAR()} of \pkg{vars}. See 'Details' and the manual of \code{VAR()}.
#' The default is \code{list(type = "const", lag.max = 6, ic = "AIC")}.
#' @param n.ahead An integer specifying the number of steps ahead for prediction.
#' 
#' @examples
#' library(HDTSA)
#' data(FamaFrench, package = "HDTSA")
#' 
#' ## Remove the market effects
#' reg <- lm(as.matrix(FamaFrench[, -c(1:2)]) ~ as.matrix(FamaFrench$MKT.RF))
#' Y_2d = reg$residuals
#' 
#' ## Rearrange Y_2d into a 3-dimensional array Y
#' Y = array(NA, dim = c(NROW(Y_2d), 10, 10))
#' for (tt in 1:NROW(Y_2d)) {
#'   for (ii in 1:10) {
#'     Y[tt, ii, ] <- Y_2d[tt, (1 + 10*(ii - 1)):(10 * ii)]
#'   }
#' }
#' 
#' res_cp <- CP_MTS(Y, lag.k = 20, method = "CP.Refined")
#' pred_cp_Y <- predict(res_cp, n.ahead = 1)[[1]]
#' 
#' @return \item{Y_pred}{A list of length \code{n.ahead}, where each element is a 
#' \eqn{p \times q} matrix representing the predicted values at
#' each time step.}
#' @seealso \code{\link{CP_MTS}}
#' @export
#' @method predict mtscp
predict.mtscp <- function(object, newdata = NULL, n.ahead = 10,
                          control_ARIMA = list(), control_VAR = list(), ...) {
  params <- list(...)
  con1 <- list(ic = "aic")
  con1[(namc <- names(control_ARIMA))] <- control_ARIMA
  
  con2 <- list(type = "const", lag.max = 6, ic = "AIC")
  con2[(namc <- names(control_VAR))] <- control_VAR
  
  Rank <- object$Rank
  A <- object$A
  B <- object$B
  method <- object$method
  if (missing(newdata) || is.null(newdata)) {
    newdata <- object$f
  }
  if (Rank$d == 1){
    newdata <- if (is.matrix(newdata)) {newdata}
    else if (is.vector(newdata) && length(dim(newdata)) == 1) {matrix(newdata, ncol = 1)}
    else {stop("newdata does not satisfy the matrix or vector condition")}
    f_pred <- Forecast(newdata, con1, con2, n_step = n.ahead)$f.forecast
    Y_pred <- F_pred <- list()
    for (ii in seq_len(n.ahead)) {
      F_pred[[ii]] <- f_pred[ii, ]
      Y_pred[[ii]] <- A %*% F_pred[[ii]] %*% t(B)
      }
    return(Y_pred)
    } 
  else {
    colnames(newdata) <- paste0("y", seq_len(ncol(newdata)))
    var_f <- do.call(vars::VAR, c(list(newdata), con2))

    pre_f <- predict(var_f, n.ahead = n.ahead)
    f_pred <- vector()
    for (jj in seq_len(ncol(newdata))) {
      f_pred <- cbind(f_pred, pre_f$fcst[[jj]][, 1])
    }
    colnames(f_pred)  <- paste0("f_pre", seq_len(ncol(newdata)))
    row.names(f_pred) <- paste(1:n.ahead, "step")
    Y_pred <- F_pred <- list()
    for (ii in seq_len(n.ahead)) {
      F_pred[[ii]] <- f_pred[ii, ]
      Y_pred[[ii]] <- A %*% diag(F_pred[[ii]]) %*% t(B)
    }
    return(Y_pred)
  }
}



Forecast <- function(f, con1, con2, n_step = 10){
  p <- ncol(f)
  if (p > 1) {
    colnames(f) <- paste0("y", seq_len(p))
    var_f <- do.call(vars::VAR, c(list(f), con2))
    # var_f <- vars::VAR(f, type = "const", lag.max = 6, ic = "HQ")
    pre_f <- predict(var_f, n.ahead = n_step)
    f_forecast <- vector()
    for (jj in seq_len(p)) {
      f_forecast <- cbind(f_forecast, pre_f$fcst[[jj]][, 1])
    }
    colnames(f_forecast)  <- paste0("f_pre", seq_len(p))
    row.names(f_forecast) <- paste(1:n_step, "step")
  }else {
    fit_f <- do.call(forecast::auto.arima, c(list(f), con1))
    # fit_f  <- forecast::auto.arima(f)
    pre_f <- predict(fit_f, n.ahead = n_step)
    f_forecast <- as.matrix(pre_f$pred)
    colnames(f_forecast)  <- paste0("f_pre", seq_len(p))
    row.names(f_forecast) <- paste(1:n_step, "step")
  }
  return(list(f.forecast = f_forecast))
}