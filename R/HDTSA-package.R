#' HDTSA: High Dimensional Time Series Analysis Tools
#'
#'An implementation for high-dimensional time series analysis methods, including 
#'factor model for vector time series proposed by Lam and Yao (2012)
#'<doi:10.1214/12-AOS970> and Chang, Guo and Yao (2015)
#'<doi:10.1016/j.jeconom.2015.03.024>, martingale difference test proposed by 
#'Chang, Jiang and Shao (2023) <doi:10.1016/j.jeconom.2022.09.001>, principal 
#'component analysis for vector time series proposed by Chang, Guo and Yao (2018)
#'<doi:10.1214/17-AOS1613>, cointegration analysis proposed by Zhang, Robinson
#'and Yao (2019) <doi:10.1080/01621459.2018.1458620>, unit root test proposed by
#'Chang, Cheng and Yao (2022) <doi:10.1093/biomet/asab034>, white noise test
#'proposed by Chang, Yao and Zhou (2017) <doi:10.1093/biomet/asw066> and Chang
#'et al. (2026+), CP-decomposition for matrix time series proposed by Chang et al.
#'(2023) <doi:10.1093/jrsssb/qkac011> and Chang et al. (2026)
#'<doi:10.48550/arXiv.2410.05634>, and statistical inference for spectral density
#'matrix proposed by Chang et al. (2025) <doi:10.1080/01621459.2025.2468013>.
#'
#' @name HDTSA-package
#' @useDynLib HDTSA
#' @importFrom Rcpp evalCpp
#' @aliases HDTSA-package HDTSA
#' @author Jinyuan Chang, Jing He, Chen lin, Qiwei Yao
#' @keywords internal
"_PACKAGE"