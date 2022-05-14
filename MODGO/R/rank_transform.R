#' Rank-based inverse normal transformation
#' 
#' Applies a rank-based inverse normal transformation to a vector.
#' 
#' @param x a numeric vector 
#' @param ties_method Method on how to deal with equal values 
#' during rank transformation. Acceptable input:"max","average","min". This
#' parameter is passed to the parameter \code{ties.method} of 
#' \code{\link[base]{rank}}.
#' @return A numeric vector with the results of applying a rank-based inverse
#' normal transformation to \code{x}.
#' @author Andreas Ziegler, Francisco M. Ojeda, George Koliopanos
#' @references Beasley, T Mark, Stephen Erickson, and David B Allison (2009), 
#' ``Rank-Based Inverse Normal Transformations Are Increasingly Used, but Are
#'   They Merited?,'' Behavior Genetics 39 (5), 580–95,
#'   \url{https://doi.org/10.1007/s10519-009-9281-0}.
#' @examples
#' 
#' data("Cleveland",package="modgo")
#' test_rank <- rbi_normal_transform(Cleveland[,1])
#' 
#' @keywords Normal rank transformation
#' @export



rbi_normal_transform <- function (x, ties_method= c("max", "min", "average")) {
  
  ties_method <- match.arg(ties_method)
  
  n <- sum(!is.na(x))
  ranks <- base::rank(x, na.last = "keep", ties.method = ties_method)
  qnorm(ranks/(n + 1))
}