#' Inverse of rank-based inverse normal transformation
#' 
#' Transforms a vector \code{x} using the inverse of a rank-based inverse normal 
#' transformation associated with a given vector \code{x_original} in a 
#' sample survey. This inverse
#' is defined as \eqn{F_n^{-1}\Phi(x)}, where \eqn{F_n^{-1}} is the inverse 
#' empirical cumulative distribution function or quantile function of 
#' \code{x_original} and \eqn{\Phi} is the cumulative distribution function of 
#' a standard normal random variable.
#' 
#' @param x a numeric vector to which the inverse of a rank-based inverse normal 
#' transformation associated with \code{x_original} will be applied.
#' @param survey design object from package `survey`.
#' @param x_original name of numeric variable in in design object.  
#' @return A numeric vector.
#' @author Andreas Ziegler, Francisco M. Ojeda, George Koliopanos
#' 
#' @examples
#' 
#' library(survey)
#' data("api", package="survey")
#' # Example x_original continuous
#' dstrat <- svydesign(id = ~ 1, strata = ~ stype, weights = ~ pw, 
#'                     data = apistrat, fpc = ~ fpc)
#' dstrat <- rbi_normal_transform_svy(x = "avg.ed", design = dstrat, 
#'                                    x_transf = "avg.ed_transf")
#'                                    
#' x_transf <-  dstrat[["variables"]][["avg.ed_transf"]]  
#' x_inv_transf <- rbi_normal_transform_inv_svy(x = x_transf, design = dstrat,
#'                                              x_original = "avg.ed")
#' x_original <- dstrat[["variables"]][["avg.ed"]]                                              
#' plot(x_inv_transf, x_original)
#' 
#' # Example x_original binary
#' dstrat <- update(dstrat, yr.rnd.bin = ifelse(yr.rnd == 'Yes', 1, 0))
#' dstrat <- rbi_normal_transform_svy(x = "yr.rnd.bin", design = dstrat, 
#'                                    x_transf = "yr.rnd.bin_transf")
#' x_transf <-  dstrat[["variables"]][["yr.rnd.bin_transf"]]  
#' x_inv_transf <- rbi_normal_transform_inv_svy(x = x_transf, design = dstrat,
#'                                              x_original = "yr.rnd.bin")
#' x_original <- dstrat[["variables"]][["yr.rnd.bin"]]                                              
#' table(x_inv_transf, x_original, useNA = "ifany") 
#'                                  
#' @keywords Inverse transformation
#' @export

rbi_normal_transform_inv_svy <- function (x, design, x_original) {
  
  if (!requireNamespace("survey", quietly = TRUE)) {
    stop(
      "Package \"survey\" must be installed to use this function.",
      call. = FALSE
    )
  }    
  
  out <- rep(NA, length(x))
  x_unique <- na.omit(unique(x))
  
  if (length(x_unique) == 0) return(out)
  
  the_formula <- as.formula(paste("~ ", x_original))
  for (xval in x_unique) {
    aux <- survey::svyquantile(the_formula, design = design, 
                               quantile = pnorm(xval), 
                               qrule = "hf1", ci = FALSE, se = FALSE)
    aux <- aux[[x_original]][1]
    pos <- which(x == xval)
    out[pos] <- aux
  }
  
  out
}  

