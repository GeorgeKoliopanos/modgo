#' Inverse Generalized Poisson Distribution
#' 
#' Inverse transforms z values of a vector to simulated values driven by
#' the original dataset using Generalised Poisson distribution function
#' 
#' @param x a vector of z values 
#' @param x_original a data frame or matrix of the original dataset
#' @return A numeric vector.
#' @author Andreas Ziegler, Francisco M. Ojeda, George Koliopanos
#' 
#' @examples
#' data("Cleveland",package="modgo")
#' test_rank <- rbi_normal_transform(Cleveland[,1])
#' test_inv_rank <- poisson_transform_inv(x = test_rank,
#'                   x_original = Cleveland[,1])
#' 
#' 
#' 
#' @keywords GLDEX Inverse transformation
#' @export



poisson_transform_inv <- function (x, x_original) {
  
  if (!requireNamespace("gp", quietly = TRUE)) {
    stop(
      "Package \"gp\" must be installed to use this function.",
      call. = FALSE
    )
  }
  
  gp_var <- gp::gp.mle(x_original)
  theta <- gp_var["theta"]
  lambda <- gp_var["lambda"]
  n <- length(x)
  
  mym <- theta/(1- lambda)
  myv <- sqrt(theta*(1-lambda)^(-3))
  Y <- rnorm(n)
  var_pois <- myv^2
  if (var_pois >= 10){
    X <- round(mym+myv*Y+0.5)
  }else {X <- floor(mym+myv*Y+0.5)}
  X[X<0] <- 0
  return(X)
}