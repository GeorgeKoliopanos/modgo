#' Inverse gldex transformation
#' 
#' Inverse transforms z values of a vector to simulated values driven by
#' the original dataset using Generalized Lambda and Generalized Poisson 
#' percentile functions
#' 
#' @param x a vector of z values 
#' @param x_original a data frame or matrix of the original dataset
#' @return A numeric vector.
#' @author Andreas Ziegler, Francisco M. Ojeda, George Koliopanos
#' 
#' @examples
#' data("Cleveland",package="modgo")
#' test_rank <- rbi_normal_transform(Cleveland[,1])
#' test_inv_rank <- general_transform_inv(x = test_rank,
#'                   x_original = Cleveland[,1])
#' 
#' 
#' 
#' @keywords Generalized Inverse transformation
#' @export



general_transform_inv <- function (x, x_original, lmbds) {
  
  if (!requireNamespace("GLDEX", quietly = TRUE)) {
    stop(
      "Package \"GLDEX\" must be installed to use this function.",
      call. = FALSE
    )
  }
  
  if (length(na.omit(lmbds)) == 5){
    if(lmbds[5] == 1){model <- "fmkl"} else{model <- "rs"}
    y <- pnorm(x)
    Q <- GLDEX::qgl(y, lmbds[1:4], param = model)
  }else if (length(na.omit(lmbds)) == 2){
    theta <- lmbds[1]
    lambda <- lmbds[2]
    n <- length(x)
    
    mym <- theta/(1- lambda)
    myv <- sqrt(theta*(1-lambda)^(-3))
    Y <- rnorm(n)
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
    Q <- na.omit(unlist(GLDEX::fun.simu.bimodal(lmbds[1:4],
                                                lmbds[5:8],
                                                lmbds[9],
                                                len = length(x_original),
                                                no.test = 20,
                                                param1 = model_1,
                                                param2 = model_2)))
    Q <- as.vector(Q[1:length(x_original)])
  }else {
    stop(paste0("Error with lambda creation in ",i))
    
  }
  
  return(Q)
  
}