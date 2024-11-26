#'@method print factors
#'@export
print.factors <- function(x, digits = max(5, getOption("digits") - 5), prefix = "\t", ...){
  cat("\n")
  if(!is.null(x$reg.coff.mat)){
    title <- "Factor analysis with observed regressors for vector time series"
  }
  else{
    title <- "Factor analysis for vector time series"
  }
  cat(strwrap(title, prefix = prefix), sep = "\n")
    # cat(strwrap(x$method, prefix = prefix), sep = "\n")
  cat("\n")
  out <- character()
  if(!is.null(x$factor_num))
    if(length(x$factor_num)<2){
      out <- c(out, paste("The estimated number of factors", "=",
                          format(x$factor_num[1], digits = max(1L, digits - 2L))))
    }
    else if(length(x$factor_num)==2){
      out <- c(out, paste("The estimated number of factors in the first step", "=",
                          format(x$factor_num[1], digits = max(1L, digits - 2L))))
      out <- c(out, paste("The estimated number of factors in the second step", "=",
                          format(x$factor_num[2], digits = max(1L, digits - 2L))))
      out <- c(out, paste("The estimated number of factors", "=",
                          format(sum(x$factor_num), digits = max(1L, digits - 2L))))
    }
      
  if(!is.null(x$lag.k))
    out <- c(out, paste(names(x$lag.k), "=",
                        format(x$lag.k, digits = max(1L, digits - 2L))))
  # cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  cat(strwrap(paste(out)), sep = "\n")
  invisible(x)
}

#'@export
print.hdtstest <- function(x, digits = max(5, getOption("digits") - 5), prefix = "\t", ...){
  if(!is.null(x$type))
    title <- "Testing for martingale difference hypothesis in high dimension"
  else if(!is.null(x$cri95))
    title <- "Global hypothesis testing for spectral density matrix"
  else if(!is.null(x$MultiTest))
    title <- "Multiple testing with FDR control for spectral density matrix"
  else
    title <- "Testing for white noise hypothesis in high dimension"
  cat("\n")
  cat(strwrap(title, prefix = prefix), sep = "\n")
  cat("\n")
  out <- character()
  # if(!is.null(x$statistic))
  #   out <- c(out, paste(names(x$statistic), "=",
  #                       format(x$statistic, digits = max(1L, digits - 2L))))
  # 输出统计量和p值在一起
  if(!is.null(x$statistic) && !is.null(x$p.value)) {
    for(i in seq_along(x$statistic)) {
      out <- c(out, paste(names(x$statistic)[i], "=",
                          format(x$statistic[i], digits = max(1L, digits - 2L)), 
                          ", p-value =",
                          format.pval(x$p.value[i], digits = max(1L, digits - 3L))))
    }
  }
  if(!is.null(x$MultiTest))
    out <- c(out, paste("Resuls for multiple test", ":"), x$MultiTest
             )
  if(!is.null(x$lag.k))
    out <- c(out, paste(names(x$lag.k), "=",
                        format(x$lag.k, digits = max(1L, digits - 2L))))
  if(!is.null(x$kernel))
    out <- c(out, paste(names(x$kernel), "=",
                        format(x$kernel, digits = max(1L, digits - 2L))))
  # if(!is.null(x$p.value)) {
  #   fp <- format.pval(x$p.value, digits = max(1L, digits - 3L))
  #   out <- c(out, paste("p-value",
  #                       if(all(startsWith(fp, "<"))) fp else paste("=",fp)))
  # }
  if(!is.null(x$cri95)) {
    out <- c(out, paste(names(x$cri95), "=",
                        format(x$cri95, digits = max(1L, digits - 2L))))
  }
  
  
  # cat(strwrap(paste(out, collapse = "\n")), sep = "\n")
  # if(!is.null(x$MultiTest))
  #   cat(x$MultiTest)
  cat(strwrap(paste(out)), sep = "\n")
  if(!is.null(x$type)) {
    if(is.character(x$type))
      cat("Data map :", x$type)
      # out <- c(out, paste("Data map :",
      #                     format(x$type, digits = max(1L, digits - 2L))))
    else if(is.matrix(x$type)){
      fp <- names(x$type)[1][1]
      # out <- c(out, paste("Data map :",
      #                     format(fp, digits = max(1L, digits - 2L))))
      cat("Data map:", paste("Data map :",
                            format(fp, digits = max(1L, digits - 2L))))
    }
    else if(is.function(x$type)){
      cat("Data map:", paste(format("User", digits = max(1L, digits - 2L)),
                              collapse = ""))
      # out <- c(out, paste("Data map :",
      #                     paste(format(x$type, digits = max(1L, digits - 2L)),
      #                           collapse = "")))
    }
    # else if(typeof(x$type) == "expression"){
    #   cat("Data map:", paste(deparse(x$type[[1]])),"\n")
    # }
    # else{
    #   # out <- c(out, paste("Data map :",
    #   #                     format(x$type, digits = max(1L, digits - 2L))))
    #   # out <- c(out, paste("Data map :", x$type))
    #   cat("Data map:", paste(deparse(a)),"\n")
    # }
    
  }
  
  invisible(x)
}
#'@method print tspca
#'@export
print.tspca <- function(x, digits = max(5, getOption("digits") - 5), prefix = "\t", ...){
  cat("\n")
  cat(strwrap("Principal component analysis for time series", prefix = prefix),
      sep = "\n")
  cat("\n")
  if(!is.null(x$NoGroups)){
    cat(strwrap(x$method), sep = "\n")
    text <- "The number of groups"
    cat(paste(text, "=",
                 format(x$NoGroups, digits = max(1L, digits - 2L))),
        sep = "\n")
  }
  if(!is.null(x$No_of_Members)){
    N2 <- x$No_of_Members > 1
    text <- "The numbers of members in groups containing at least two members:"
    cat(paste(text), 
              format(x$No_of_Members[N2], digits = max(1L, digits - 2L)),sep = " ")
    cat("\n")
  }
  # if(!is.null(x$Groups)){
  #   text <- "The indices of elements in each group containing at least two elements:"
  #   cat(paste(text, "\n"))
  #   print(x$Groups[,N2])
  # }
  if(length(x) <= 4) cat(strwrap(x$method[2]), sep = "\n")
  invisible(x)
}

#'@method print mtscp
#'@export
print.mtscp <- function(x, digits = max(5, getOption("digits") - 5), prefix = "\t", ...){
  cat("\n")
  cat(strwrap("Estimation of matrix CP-factor model", prefix = prefix), sep = "\n")
  cat("\n")
  if(!is.null(x$Rank)){
    cat(strwrap(paste("Method:", x$method)), sep = "\n")
    if(length(x$Rank)>1)
      text <- c("The estimated number of factors d",
                "The estimated rank of the left loading matrix d1",
                "The estimated rank of the right loading matrix d2")
    else text <- c("The estimated number of factors d")
    cat(strwrap(paste(text,"=",format(x$Rank, digits = max(1L, digits - 2L)))),
        sep = "\n")
        
  }
  # cat(strwrap("Use $A or $B or $f to access the corresponding data in the result list"), sep = "\n")
  invisible(x)
}

#'@method print coint
#'@export
print.coint <- function(x, digits = max(5, getOption("digits") - 5), prefix = "\t", ...){
  cat("\n")
  cat(strwrap("Cointegration analysis for vector time series", prefix = prefix),
      sep = "\n")
  cat("\n")
  if(is.matrix(x$coint_rank)){
    cat(strwrap(paste("Using", x$method)), sep = "\n")
    cat(strwrap("The estimated number of cointegration rank:"), sep = "\n")
    print(x$coint_rank)
    if(!is.null(x$lag.k))
      out <- paste(names(x$lag.k), "=",
                          format(x$lag.k, digits = max(1L, digits - 2L)))
    cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  }
  else{
    cat(strwrap(paste("Using", x$method)), sep = "\n")
    out <- character()
    if(!is.null(x$coint_rank))
      out <- c(out, paste(names(x$coint_rank), "=",
                          format(x$coint_rank, digits = max(1L, digits - 2L))))
    if(!is.null(x$lag.k))
      out <- c(out, paste(names(x$lag.k), "=",
                          format(x$lag.k, digits = max(1L, digits - 2L))))
    cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  }

  invisible(x)
}

#'@method print urtest
#'@export
print.urtest <- function(x, digits = max(5, getOption("digits") - 5), prefix = "\t", ...){
  title <- "Testing for unit roots based on sample autocovariances"
  cat("\n")
  cat(strwrap(title, prefix = prefix), sep = "\n")
  cat("\n")
  if(length(x$reject)>1){
    cat(strwrap("Reject the null hypothesis or not with different argument"),
        sep = "\n")
    print(t(x$reject))
  }
  else{
    out <- character()
    if(!is.null(x$reject))
      out <- c(out, paste("reject", "=",
                          format(x$reject[1,1], digits = max(1L, digits - 2L))))
    if(!is.null(x$statistic))
      out <- c(out, paste("Statistic", "=",
                          format(x$statistic, digits = max(1L, digits - 2L))))
    cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  }
  
  invisible(x)
}
