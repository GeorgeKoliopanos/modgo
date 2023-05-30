#' Inverse of rank based inverse normal transformation
#' 
#' Transforms a vector \code{x} using the inverse of rank based inverse normal 
#' transformation associated with a given vector \code{x_original}. This inverse
#' is defined as \eqn{F_n^{-1}\Phi(x)}, where \eqn{F_n^{-1}} is the inverse 
#' empirical cumulative distribution function of \code{x_original} and 
#' \eqn{\Phi} is the cumulative distribution function of a standard normal 
#' random variable.
#' 
#' @param x a numeric vector. 
#' @param x_original a numeric vector from the original dataset
#' @return A numeric vector.
#' @author Andreas Ziegler, Francisco M. Ojeda, George Koliopanos
#' 
#' @examples
#' data("Cleveland",package="modgo")
#' test_rank <- rbi_normal_transform(Cleveland[,1])
#' test_inv_rank <- rbi_normal_transform_inv(x = test_rank,
#'                                           x_original = Cleveland[,1])
#' 
#' 
#' 
#' @keywords Inverse transformation
#' @export

 

rbi_normal_transform_inv <- function (x, x_original) {
  #type = 1 in quantile gives the inverse of the empirical cdf of x_original.
  quantile(x_original, probs = pnorm(x), type = 1, na.rm = TRUE, names = FALSE)
}