#' Rank normal transformation
#' 
#' Transforms values of a vector to ranks 
#' and ranks to z values from a normal distribution
#' 
#' @param x a numeric vector 
#' @param ties_method Method on how to deal with equal values 
#' during rank transformation.Acceptable input:"max","average","min".
#' @return A numeric vector.
#' @author Andreas Ziegler, Francisco M. Ojeda, George Koliopanos
#' @examples
#' 
#' data("Cleveland",package="modgo")
#' test_rank <- rbi_normal_transform(Cleveland[,1])
#' 
#' @keywords Normal rank transformation
#' @export



rbi_normal_transform <- function (x,ties_method="max") {
  n <- sum(!is.na(x))
  ranks <- base::rank(x, na.last = "keep", ties.method = ties_method)
  qnorm(ranks/(n + 1))
}