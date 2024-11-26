#' Make predictions from a \code{"factors"} object.
#'
#' This function makes predictions from a \code{"factors"} object.
#' 
#' @details
#' Forecasting for \eqn{{\bf y}_t} can be implemented in two steps:
#' 
#' \emph{Step 1}. Get the \eqn{h}-step-ahead forecast of the \eqn{\hat{r} \times 1} 
#' time series \eqn{\hat{\bf x}_t} [See \code{\link{Factors}}] by using a VAR model 
#' (if \eqn{\hat{r} > 1}) with the order determined by AIC or an ARIMA model (if \eqn{\hat{r} = 1}) with the order 
#' determined by AIC. The orders of VAR and ARIMA models can also be specified by 
#' users through the arguments \code{control_VAR} and \code{control_ARIMA}, respectively.
#' 
#' \emph{Step 2}. The forecasted value for \eqn{{\bf y}_t} is obtained by
#'  \eqn{\hat{\bf y}_{n+h}= \hat{\bf A}\hat{\bf x}_{n+h}}.
#' 
#' 
#'
#' @param object An object of class \code{"factors"} constructed by \code{\link{Factors}}.
#' @param newdata A new data matrix to predict from (optional).
#' @param ... Currently not used.
#' @param control_ARIMA A list of arguments passed to the function
#' \code{auto.arima()} of \pkg{forecast} [See the manual of \code{auto.arima()}].
#' See 'Details'. The default is \code{list(ic = "aic")}.
#' @param control_VAR A list of arguments passed to the function
#' \code{VAR()} of \pkg{vars} [See the manual of \code{VAR()}]. See 'Details'.
#' The default is \code{list(type = "const", lag.max = 6, ic = "AIC")}.
#' @param n.ahead An integer specifying the number of steps ahead to predict.
#' 
#' @return \item{ts_pred}{A matrix of predicted values}
#' @seealso \code{\link{Factors}}
#' @importFrom vars VAR
#' @importFrom forecast auto.arima
#' @importFrom stats predict
#' @export
predict.factors <- function(object, newdata = NULL, n.ahead = 10,
                            control_ARIMA = list(), control_VAR = list(), ...) {
  params <- list(...)
  con1 <- list(ic = "aic")
  con1[(namc <- names(control_ARIMA))] <- control_VAR
  
  con2 <- list(type = "const", lag.max = 6, ic = "AIC")
  con2[(namc <- names(control_VAR))] <- control_VAR
  
  if (!is.null(object$loading.mat)) {
    loading <- object$loading.mat
  }
  if (missing(newdata) || is.null(newdata)) {
    ts <- object$X
  }
  else {ts <- newdata %*% as.matrix(loading)}
  ts_pred <- forecast(ts, con1, con2, n.ahead)$f.forecast %*% t(loading)
  ts_pred
}

#' Make predictions from a \code{"tspca"} object.
#'
#' This function makes predictions from a \code{"tspca"} object.
#' 
#' @details
#' Having obtained \eqn{\hat{\bf x}_t} using the \code{\link{PCA_TS}} function, which is
#' segmented into \eqn{q} uncorrelated subseries
#' \eqn{\hat{\bf x}_t^{(1)},\ldots,\hat{\bf x}_t^{(q)}}, the forecasting for \eqn{{\bf y}_t}
#' can be performed in two steps:
#' 
#' \emph{Step 1}. Get the \eqn{h}-step-ahead forecast \eqn{\hat{\bf x}_{n+h}^{(j)}} \eqn{(j=1,\ldots,q)}
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
#' @param newdata A new data matrix to predict from (optional).
#' @param ... Currently not used.
#' @param control_ARIMA A list of arguments passed to the function
#' \code{auto.arima()} of \pkg{forecast} [See the manual of \code{auto.arima()}].
#' See 'Details'. The default is \code{list(max.d = 0, max.q = 0, ic = "aic")}.
#' @param control_VAR A list of arguments passed to the function
#' \code{VAR()} of \pkg{vars} [See the manual of \code{VAR()}]. See 'Details'.
#' The default is \code{list(type = "const", lag.max = 6, ic = "AIC")}.
#' @param n.ahead An integer specifying the number of steps ahead to predict.
#' 
#' @seealso \code{\link{PCA_TS}}
#' 
#' @return \item{Y_pred}{A matrix of predicted values.}
#' @export
predict.tspca <- function(object, newdata = NULL, n.ahead = 10,
                          control_ARIMA = list(), control_VAR = list(), ...) {
  params <- list(...)
  con1 <- list(max.d = 0, max.q = 0, ic = "aic")
  con1[(namc <- names(control_ARIMA))] <- control_VAR
  
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
    ts_i_pred <- forecast(ts_i, con1, con2, n.ahead)$f.forecast
    ts_pred[1:n.ahead, group[, ii]] <- ts_i_pred
  }
  Y_pred <- ts_pred %*% t(solve(B))
  Y_pred
}

#' Make predictions from a \code{"mtscp"} object.
#'
#' This function makes predictions from a \code{"mtscp"} object.
#' 
#' @details
#' Forecasting for \eqn{{\bf y}_t} can be implemented in two steps:
#' 
#' \emph{Step 1}. Get the \eqn{h}-step ahead forecast of the \eqn{\hat{d} \times 1} 
#' time series \eqn{\hat{\bf x}_t} [See \code{\link{CP_MTS}}] by using a VAR model 
#' (if \eqn{\hat{d} > 1}) or an ARIMA model
#' (if \eqn{\hat{d} = 1}). The orders of VAR and ARIMA models are determined by
#' AIC by default. Otherwise, they can also be specified by users through the
#' arguments \code{control_VAR} and \code{control_ARIMA}, respectively.
#' 
#' \emph{Step 2}. The forecasted value for \eqn{{\bf y}_t} is obtained by
#'  \eqn{\hat{\bf y}_{n+h}= \hat{\bf A} {\rm diag}(\hat{\bf x}_{n+h}) \hat{\bf B}'} .
#'
#' @param object An object of class \code{"mtscp"} constructed by \code{\link{CP_MTS}}.
#' @param newdata A new data matrix to predict from (optional).
#' @param ... Currently not used.
#' @param control_ARIMA A list of arguments passed to the function
#' \code{auto.arima()} of \pkg{forecast} [See the manual of \code{auto.arima()}].
#' See 'Details'. The default is \code{list(ic = "aic")}.
#' @param control_VAR A list of arguments passed to the function
#' \code{VAR()} of \pkg{vars} [See the manual of \code{VAR()}]. See 'Details'.
#' The default is \code{list(type = "const", lag.max = 6, ic = "AIC")}.
#' @param n.ahead An integer specifying the number of steps ahead to predict.
#' 
#' 
#' @return \code{Y_pred} A list of length \code{n.ahead}, where each element is a 
#' \eqn{p \times q} matrix representing the predicted values for 
#' each time steps.
#' @seealso \code{\link{CP_MTS}}
#' @export
predict.mtscp <- function(object, newdata = NULL, n.ahead = 10,
                          control_ARIMA = list(), control_VAR = list(), ...) {
  params <- list(...)
  con1 <- list(ic = "aic")
  con1[(namc <- names(control_ARIMA))] <- control_VAR
  
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
    # ar_f <- ar(newdata)
    # pre_f <- predict(ar_f, n.ahead = n.ahead)
    # f_pred <- as.matrix(pre_f$pred)
    # colnames(f_pred) <- paste0("f_pre", seq_len(ncol(newdata)))
    # row.names(f_pred) <- paste(1:n.ahead, "step")
    f_pred <- forecast(newdata, con1, con2, n_step = n.ahead)$f.forecast
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
    # var_f <- try(vars::VAR(newdata, type = "const", lag.max = 6, ic = "HQ"),
    #              silent = TRUE)

    # if (inherits(var_f, "try-error")) {
    #   f_pred <- vector()
    # 
    #   for (jj in seq_len(Rank$d)){
    #     fit_f <- do.call(forecast::auto.arima, c(list(newdata[, jj]), con1))
    #     # fit_f <- forecast::auto.arima(newdata[, jj])
    #     # fit_f <- ar(newdata[, jj])
    #     pre_f <- predict(fit_f, n.ahead = n.ahead)
    #     f_pred_i <- as.matrix(pre_f$pred)
    #     f_pred <- cbind(f_pred, f_pred_i)
    #   }
    # 
    #   colnames(f_pred) <- paste0("f_pre", seq_len(ncol(newdata)))
    #   row.names(f_pred) <- paste(1:n.ahead, "step")
    # }

    pre_f <- predict(var_f, n.ahead = n.ahead)
    f_pred <- vector()
    for (jj in seq_len(ncol(newdata))) {
      f_pred <- cbind(f_pred, pre_f$fcst[[jj]][, 1])
    }
    colnames(f_pred)  <- paste0("f_pre", seq_len(ncol(newdata)))
    row.names(f_pred) <- paste(1:n.ahead, "step")
    # if (c(is.na(f_pred))[1]) {
    #   f_pred <- vector()
    # 
    #   for (jj in seq_len(Rank$d)){
    #     fit_f <- do.call(forecast::auto.arima, c(list(newdata[, jj]), con1))
    #     pre_f <- predict(fit_f, n.ahead = n.ahead)
    #     f_pred_i <- as.matrix(pre_f$pred)
    #     f_pred <- cbind(f_pred, f_pred_i)
    #   }
    # 
    #   colnames(f_pred) <- paste0("f_pre", seq_len(ncol(newdata)))
    #   row.names(f_pred) <- paste(1:n.ahead, "step")
    #   }
    Y_pred <- F_pred <- list()
    for (ii in seq_len(n.ahead)) {
      F_pred[[ii]] <- f_pred[ii, ]
      Y_pred[[ii]] <- A %*% diag(F_pred[[ii]]) %*% t(B)
    }
    return(Y_pred)
  }
}



forecast <- function(f, con1, con2, n_step = 10){
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