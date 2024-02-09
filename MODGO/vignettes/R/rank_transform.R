#' Rank based inverse normal transformation
#' 
#' Applies the rank based inverse normal transformation to numeric vector.
#' 
#' The rank based inverse normal transformation (Beasley et al. (2009)), 
#' transforms values of a vector to ranks and then applies the quantile 
#' function of the standard normal distribution.
#' 
#' @param x a numeric vector 
#' @param ties_method Method on how to deal with equal values 
#' during rank transformation.Acceptable input:"max","average","min". This
#' parameter is passed to the parameter \code{ties.method} of 
#' \code{\link[base]{rank}}.
#' @return A numeric vector.
#' @author Andreas Ziegler, Francisco M. Ojeda, George Koliopanos
#' @references Beasley, T.M. and Erickson S. and Allison D.B. (2009), 
#' ``Rank-based inverse normal transformations are increasingly used, 
#' but are they merited?,'' \emph{Behavior genetics} \strong{39}, 580-595.
#' @examples
#' 
#' data("Cleveland",package="modgo")
#' test_rank <- rbi_normal_transform(Cleveland[,1])
#' 
#' @keywords Normal rank transformation
#' @export
#' @import stats


rbi_normal_transform <- function (x, ties_method= c("max", "min", "average")) {
  
  ties_method <- match.arg(ties_method)
  
  n <- sum(!is.na(x))
  ranks <- base::rank(x, na.last = "keep", ties.method = ties_method)
  qnorm(ranks/(n + 1))
}