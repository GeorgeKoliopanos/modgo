#' Inverse gldex transformation
#' 
#' Inverse transforms z values of a vector to simulated values driven by
#' the original dataset using Generalized Lambda and Generalized Poisson 
#' percentile functions.
#' 
#' @param x A vector of z values.
#' @param data A data frame with original variables.
#' @param n_samples Number of samples you need to produce.
#' @param lmbds A vector with generalized lambdas values
#' @return A numeric vector.
#' @author Andreas Ziegler, Francisco M. Ojeda, George Koliopanos
#' 
#' @examples
#' data("Cleveland",package="modgo")
#' test_rank <- rbi_normal_transform(Cleveland[,1])
#' test_generalized_lmbds <- generalizedMatrix(Cleveland, 
#'                   bin_variables = c("Sex", "HighFastBloodSugar", "CAD"))
#' test_inv_rank <- general_transform_inv(x = test_rank,
#'                   data = Cleveland[,1],
#'                   n_samples = 100,
#'                   lmbds = test_generalized_lmbds[,1])
#' 
#' 
#' 
#' @keywords Generalized Inverse transformation
#' @export
#' @import GLDEX
#' @import gp
#' @import stats 



general_transform_inv <- function (x,
                                   data = NULL,
                                   n_samples,
                                   lmbds) {
  
  if (length(na.omit(lmbds)) == 5){
    if(lmbds[5] == 1){model <- "fmkl"} else{model <- "rs"}
    y <- pnorm(x)
    Q <- GLDEX::qgl(y, lmbds[1:4], param = model)
  }else if (length(na.omit(lmbds)) == 1){
    Q <- rbi_normal_transform_inv(x,
                                  rbinom(n = length(x),
                                         1,
                                         prob = lmbds[1]))
  }else if (length(na.omit(lmbds)) == 2){
    theta <- lmbds[1]
    lambda <- lmbds[2]

    mym <- theta/(1- lambda)
    myv <- sqrt(theta*(1-lambda)^(-3))
    Y <- x
    var_pois <- myv^2
    if (var_pois >= 10){
      Q <- round(mym+myv*Y+0.5)
    }else {
      Q <- floor(mym+myv*Y+0.5)
      }
    Q[Q<0] <- 0
  
  }else if (length(na.omit(lmbds)) == 11){
    if(lmbds[10] == 1){model_1 <- "fmkl"} else{model_1 <- "rs"}
    if(lmbds[11] == 1){model_2 <- "fmkl"} else{model_2 <- "rs"}
    
    first_distr_num <- ceiling(n_samples*lmbds[9])
    y <- pnorm(x)
    Q_1 <- GLDEX::qgl(y[1:first_distr_num],
                      lmbds[1:4],
                      param = model_1)
    Q_2 <- GLDEX::qgl(y[(first_distr_num + 1):n_samples],
                      lmbds[5:8],
                      param = model_2)
    Q <- as.vector(c(Q_1,Q_2))
  }else if (length(na.omit(lmbds)) == 0){
    Q <- rbi_normal_transform_inv(x, data)
    
  }
  
  return(Q)
  
}