#' Inverse gldex transformation
#' 
#' Inverse transforms z values of a vector to simulated values driven by
#' the original dataset using Generalised Lambda percentile function
#' 
#' @param x a vector of z values 
#' @param x_original a data frame or matrix of the original dataset
#' @return A numeric vector.
#' @author Andreas Ziegler, Francisco M. Ojeda, George Koliopanos
#' 
#' @examples
#' data("Cleveland",package="modgo")
#' test_rank <- rbi_normal_transform(Cleveland[,1])
#' test_inv_rank <- gldex_transform_inv(x = test_rank,
#'                   x_original = Cleveland[,1])
#' 
#' 
#' 
#' @keywords GLDEX Inverse transformation
#' @export



gldex_transform_inv <- function (x, x_original) {
  
  if (!requireNamespace("GLDEX", quietly = TRUE)) {
    stop(
      "Package \"GLDEX\" must be installed to use this function.",
      call. = FALSE
    )
  }
  
  lmbds <- GLDEX::fun.RMFMKL.ml(x_original)
  y <- pnorm(x)
  Q <- GLDEX::qgl(y, lmbds, param = "fmkl")
}