#'@method print factors
#'@export
print.factors <- function(x, digits = max(5, getOption("digits") - 5), prefix = "\t", ...){
  cat("\n")
  cat(strwrap(x$method, prefix = prefix), sep = "\n")
  cat("\n")
  out <- character()
  if(!is.null(x$factor_num))
    out <- c(out, paste(names(x$factor_num), "=",
                        format(x$factor_num, digits = max(1L, digits - 2L))))
  if(!is.null(x$lag.k))
    out <- c(out, paste(names(x$lag.k), "=",
                        format(x$lag.k, digits = max(1L, digits - 2L))))
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  invisible(x)
}

#'@export
print.hdtstest <- function(x, digits = max(5, getOption("digits") - 5), prefix = "\t", ...){
  cat("\n")
  cat(strwrap(x$method, prefix = prefix), sep = "\n")
  cat("\n")
  out <- character()
  if(!is.null(x$statistic))
    out <- c(out, paste(names(x$statistic), "=",
                        format(x$statistic, digits = max(1L, digits - 2L))))
  if(!is.null(x$MultiTest))
    out <- c(out, paste("Resuls for multiple test", ":")
             )
  if(!is.null(x$lag.k))
    out <- c(out, paste(names(x$lag.k), "=",
                        format(x$lag.k, digits = max(1L, digits - 2L))))
  if(!is.null(x$kernel))
    out <- c(out, paste(names(x$kernel), "=",
                        format(x$kernel, digits = max(1L, digits - 2L))))
  if(!is.null(x$p.value)) {
    fp <- format.pval(x$p.value, digits = max(1L, digits - 3L))
    out <- c(out, paste("p-value",
                        if(startsWith(fp, "<")) fp else paste("=",fp)))
  }
  if(!is.null(x$cri95)) {
    out <- c(out, paste(names(x$cri95), "=",
                        format(x$cri95, digits = max(1L, digits - 2L))))
  }
  if(!is.null(x$type)) {
    if(is.character(x$type))
      out <- c(out, paste("Data map :",
                          format(x$type, digits = max(1L, digits - 2L))))
    else if(is.matrix(x$type)){
      fp <- names(x$type)[1][1]
      out <- c(out, paste("Data map :",
                         format(fp, digits = max(1L, digits - 2L))))
    }
    else if(is.function(x$type)){
      out <- c(out, paste("Data map :",
                          paste(format(x$type, digits = max(1L, digits - 2L)),
                                collapse = "")))
    }
    else{
      out <- c(out, paste("Data map :",
                          format(x$type, digits = max(1L, digits - 2L))))
    }
    
  }
  
  
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  if(!is.null(x$MultiTest))
    cat(x$MultiTest)
  invisible(x)
}
#'@method print tspca
#'@export
print.tspca <- function(x, digits = max(5, getOption("digits") - 5), prefix = "\t", ...){
  cat("\n")
  cat(strwrap(x$method[1], prefix = prefix), sep = "\n")
  cat("\n")
  if(!is.null(x$NoGroups)){
    cat(strwrap(x$method[2]), sep = "\n")
    text <- "The number of groups"
    cat(paste(text, "=",
                 format(x$NoGroups, digits = max(1L, digits - 2L))),
        sep = "\n")
  }
  if(!is.null(x$No_of_Members)){
    N2 <- x$No_of_Members > 1
    text <- "The number of members in each of groups with at least two members"
    cat(paste(text, ":"), 
              format(x$No_of_Members[N2], digits = max(1L, digits - 2L)),sep = " ")
    cat("\n")
  }
  if(!is.null(x$Groups)){
    text <- "The indices of components in each of groups with at least two members"
    cat(paste(text, ":\n"))
    print(x$Groups[,N2])
  }
  if(length(x) <= 4) cat(strwrap(x$method[2]), sep = "\n")
  invisible(x)
}

#'@method print mtscp
#'@export
print.mtscp <- function(x, digits = max(5, getOption("digits") - 5), prefix = "\t", ...){
  cat("\n")
  cat(strwrap(x$method[1], prefix = prefix), sep = "\n")
  cat("\n")
  if(!is.null(x$Rank)){
    cat(strwrap(x$method[2]), sep = "\n")
    cat(strwrap(paste(names(x$Rank),"=",format(x$Rank, digits = max(1L, digits - 2L)))),
        sep = "\n")
        
  }
  cat(strwrap("Use $A or $B or $f to access the corresponding data in the result list"), sep = "\n")
  invisible(x)
}

#'@method print coint
#'@export
print.coint <- function(x, digits = max(5, getOption("digits") - 5), prefix = "\t", ...){
  cat("\n")
  cat(strwrap(x$method[1], prefix = prefix), sep = "\n")
  cat("\n")
  if(is.matrix(x$coint_rank)){
    cat(strwrap(x$method[2]), sep = "\n")
    cat(strwrap("The estimated number of cointegration rank:"), sep = "\n")
    print(x$coint_rank)
    if(!is.null(x$lag.k))
      out <- paste(names(x$lag.k), "=",
                          format(x$lag.k, digits = max(1L, digits - 2L)))
    cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  }
  else{
    cat(strwrap(x$method[2]), sep = "\n")
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
  cat("\n")
  cat(strwrap(x$method, prefix = prefix), sep = "\n")
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
