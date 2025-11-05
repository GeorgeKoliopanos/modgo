#' Rank-based inverse normal transformation for sample surveys
#' 
#' Applies a weighted rank-based inverse normal transformation to a numeric 
#' vector.
#' 
#' For a simple random sample the rank-based inverse normal transformation 
#' (Beasley et al. (2009)), transforms values of a vector x of length n to 
#' ranks/(n + 1) 
#' and then applies the quantile function of the standard normal distribution. 
#' An analog for sample surveys can be obtained based on the observation that 
#' in a simple random sample the ranks of x can be obtained as n times the 
#' cumulative empirical distribution of x evaluated at x. 
#' For a complex survey the Hajek 
#' estimator of the cumulative distribution function is used multiplied by
#' \eqn{\frac{\hat{N}}{\hat{N} + 1}}, where \eqn{\hat{N}} is the 
#' Horvitz-Thompson estimator of the population size.
#' 
#' @param x character string giving name of numeric variable in design object  
#' @param survey design object from survey package.
#' @param x_transf character string giving name of variable in design object
#' that will contain the rank transformed x.  
#' @return The survey design object with the added rank-transformed variable.
#' @author Andreas Ziegler, Francisco M. Ojeda, George Koliopanos
#' @references Beasley, T.M. and Erickson S. and Allison D.B. (2009). 
#' Rank-based inverse normal transformations are increasingly used, 
#' but are they merited? \emph{Behavior genetics}, 39, 580-595.
#' @examples
#' 
#' library(survey)
#' data("api", package="survey")
#' dstrat <- svydesign(id = ~ 1, strata = ~ stype, weights = ~ pw, 
#'                     data = apistrat, fpc = ~ fpc)
#' dstrat <- rbi_normal_transform_svy(x = "avg.ed", design = dstrat, 
#'                                    x_transf = "avg.ed_transf")
#' svyhist(~ avg.ed, dstrat) 
#' svyhist(~ avg.ed_transf, dstrat)                                    
#' @keywords Normal rank transformation
#' @export


rbi_normal_transform_svy <- function (x, design, 
                                      x_transf = paste0(x, "_rbnint")) {
  
  
  if (!requireNamespace("survey", quietly = TRUE)) {
    stop(
      "Package \"survey\" must be installed to use this function.",
      call. = FALSE
    )
  }  
  
  # empirical cumulative distribution function
  wgt_cdf <- survey::svycdf(as.formula(paste("~", x)), design = design)
  
  w <- weights(design, type = "sampling")
  n_total <- sum(w)
  ranks_wgt <- wgt_cdf[[1]](design[["variables"]][[x]])
  ranks_wgt <- ranks_wgt * n_total / (n_total + 1)
  rbi_result <- qnorm(ranks_wgt)
  
  design[["variables"]][[x_transf]] <- rbi_result
  design
  
}

