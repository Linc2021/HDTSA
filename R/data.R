#' U.S. Industrial Production indices
#'
#' @description  The dataset consists of 7 monthly U.S. Industrial
#' Production indices, namely \emph{the total index}, \emph{nonindustrial supplies},
#' \emph{final products}, \emph{manufacturing}, \emph{materials}, \emph{mining},
#' and \emph{utilities}, from January 1947 to December 2023 published by the U.S.
#' Federal Reserve.
#' 
#' @name IPindices
#' @docType data
#' @keywords data
#' @usage data(IPindices)
#' @source \url{https://fred.stlouisfed.org/release/tables?rid=13&eid=49670}
#'
#' @rdname US-Industrial-Production-indices
#' @format A data frame with 924 rows and 8 variables:
#' \describe{
#'   \item{DATE}{The observation date}
#'   \item{INDPRO}{The total index}
#'   \item{IPB54000S}{Nonindustrial supplies}
#'   \item{IPFINAL}{Final products}
#'   \item{IPMANSICS}{Manufacturing}
#'   \item{IPMAT}{Materials}
#'   \item{IPMINE}{Mining}
#'   \item{IPUTIL}{Utilities}
#' }
"IPindices"

#' Fama-French 10*10 return series
#'
#' @description  The portfolios are constructed by the intersections of 10 levels
#' of size, denoted by \eqn{{\rm S}_{1}, \ldots, {\rm S}_{10}}, and 10 levels of the book
#' equity to market equity ratio (BE), denoted by \eqn{{\rm BE}_1, \ldots,{\rm BE}_{10}}.
#' The dataset consists of monthly returns from January 1964 to
#' December 2021, which contains 69600 observations for 696 total months.
#' 
#' @name FamaFrench
#' @docType data
#' @keywords data
#' @usage data(FamaFrench)
#' @source \url{http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html}
#'
#' @rdname Fama-French-data
#' @format A data frame with 696 rows and 102 columns. The first column
#' represents the month, and the second column named
#' \code{MKT.RF} represents the monthly market returns. The rest of the columns
#' represent the return series for different sizes and BE-ratios. 
"FamaFrench"

#' The national QWI hires data
#'
#' @description  The data on new hires at a national level are obtained from the
#' Quarterly Workforce Indicators (QWI) of the Longitudinal Employer-Household
#' Dynamics program at the U.S. Census Bureau (Abowd et al., 2009). The national
#' QWI hires data covers a variable number of years, with some states providing
#' time series going back to 1990 (e.g., Washington), and others
#' (e.g., Massachusetts) only commencing at 2010. For each of 51 states
#' (excluding D.C. but including Puerto Rico) there is a new hires time series
#' for each county. Additional description of the data, along with its relevancy
#' to labor economics, can be found in Hyatt and McElroy (2019).
#' 
#' @name QWIdata
#' @docType data
#' @keywords data
#' @usage data(QWIdata)
#' @source \url{https://qwiexplorer.ces.census.gov/}
#' 
#' \url{https://ledextract.ces.census.gov/qwi/all}
#' @references 
#' Abowd, J. M., Stephens, B. E., Vilhuber, L., Andersson, F., McKinney, K. L., Roemer,
#' M., and Woodcock, S. (2009). The LEHD infrastructure files and the creation of the
#' quarterly workforce indicators. In \emph{Producer dynamics: New evidence from micro data},
#' pages 149--230. University of Chicago Press. \doi{10.7208/chicago/9780226172576.003.0006}.
#' 
#' Hyatt, H. R. and McElroy, T. S. (2019). Labor reallocation, employment,
#'  and earnings: Vector autoregression evidence. \emph{Labour}, \strong{33}, 463--487.
#'  \doi{doi:10.1111/labr.12153}
#'
#' @rdname QWIdata
#' @format A list with 51 elements. Every element contains a multivariate time series.
"QWIdata"

