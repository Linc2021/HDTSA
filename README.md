An introduction to `HDTSA`
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

# HDTSA

<!-- badges: start -->
<!-- badges: end -->

An implementation for high-dimensional time series analysis
methods,including: factor model for vector time series proposed by Lam
and Yao (2012)
[\<doi:10.1214/12-AOS970\>](https://doi.org/10.1214/12-AOS970) and
Chang, Guo and Yao (2015)
[\<doi:10.1016/j.jeconom.2015.03.024\>](https://doi.org/10.1016/j.jeconom.2015.03.024),
martingale difference test proposed by Chang, Jiang and Shao (2022)
[\<doi:10.1016/j.jeconom.2022.09.001\>](https://doi.org/10.1016/j.jeconom.2022.09.001),
principal component analysis for vector time series proposed by Chang,
Guo and Yao (2018)
[\<doi:10.1214/17-AOS1613\>](https://doi.org/10.1214/17-AOS1613),
cointegration analysis proposed by Zhang, Robinson and Yao (2019)
[\<doi:10.1080/01621459.2018.1458620\>](https://doi.org/10.1080/01621459.2018.1458620),
unit root test proposed by Chang, Cheng and Yao (2021)
[\<doi:10.1093/biomet/asab034\>](https://doi.org/10.1093/biomet/asab034),
white noise test proposed by Chang, Yao and Zhou (2017)
[\<doi:10.1093/biomet/asw066\>](https://doi.org/10.1093/biomet/asw066),
CP-decomposition for high-dimensional matrix time series proposed by
Chang, He, Yang and Yao (2023)
[\<doi:10.1093/jrsssb/qkac011\>](https://doi.org/10.1093/jrsssb/qkac011)
and Chang, Du, Huang and Yao (2024)
[\<doi:10.48550/arXiv.2410.05634\>](https://doi.org/10.48550/arXiv.2410.05634),
and Statistical inference for high-dimensional spectral density matrix
porposed by Chang, Jiang, McElroy and Shao (2023)
[\<doi:10.48550/arXiv.2212.13686\>](https://doi.org/10.48550/arXiv.2212.13686).

## Installation

You can install the released version of `HDTSA` from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("HDTSA")
```

Or try the development version on GitHub:

``` r
# install.packages("devtools")
devtools::install_github("Linc2021/HDTSA")
```

## Example

This is a basic example which shows you how to solve a unit root test
problem :

``` r
library(HDTSA)
N=100
Y=arima.sim(list(ar=c(0.9)), n = 2*N, sd=sqrt(1))
con_vec=c(0.45,0.55,0.65)
lagk.vec=c(0,1,2)
UR_test(Y,lagk.vec=lagk.vec, con_vec=con_vec,alpha=0.05)
#> 
#>  Testing for unit roots based on sample autocovariances
#> 
#> Reject the null hypothesis or not with different argument
#>          time_lag=0 time_lag=1 time_lag=2
#> con=0.45          0          0          0
#> con=0.55          0          0          0
#> con=0.65          0          0          0
```

``` r
UR_test(Y,alpha=0.05)
#> 
#>  Testing for unit roots based on sample autocovariances
#> 
#> Reject the null hypothesis or not with different argument
#>          time_lag=0 time_lag=1 time_lag=2 time_lag=3 time_lag=4
#> con=0.55          0          0          0          0          0
```

Here, we have provided just one example. You can use functions within
the package `HDTSA` to solve other problems. For details, please refer
to

``` r
help("HDTSA")
```

## Bug report

Please send an email to Chen Lin(<linchen@smail.swufe.edu.cn>).
