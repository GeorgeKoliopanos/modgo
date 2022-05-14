#' Correlation transformation
#' 
#' Transforming covariance matrix using biserial, tetrachoric and
#' polyserial functions 
#' 
#' This function is used by internally by \code{\link[modgo]{modgo}}.
#' 
#' @param data a data frame containing the original data.
#' @param data_z a data frame where the continuous variables in \code{data} 
#' have transformed using the rank-based inverse normal transformation. Other
#' variables are unchanged. 
#' @param Sigma A square matrix of size \code{length(variables)} containing
#' the covariances between the rank-based inverse normal transformed continuous
#' variables, and ones in the diagonal. 
#' @param variables a character vector listing the variables in the covariance
#' matrix \code{Sigma}
#' @param bin_variables  a character vector listing the binary variables in 
#' \code{variables}.
#' @param categ_variables a character vector listing the ordinal categorical in 
#' \code{variables}.
#' @return The matrix \code{Sigma} with all entries corresponding to binary
#' and categorical ordinal variables filled using the corresponding correlation
#' function (biserial, tetrachoric, polyserial).
#' @author Francisco M. Ojeda, George Koliopanos
#' @keywords Normal rank transformation
#' @examples
#'
#'@importFrom psych biserial
#'@importFrom psych polychoric
#'@importFrom psych tetrachoric
#'@importFrom psych polyserial

Sigma_transformation <- function(data,data_z,Sigma,variables,
                                 bin_variables=c(),categ_variables=c()) {

  cont_variables <- setdiff(variables,c(bin_variables,categ_variables))
  #Binary correlations
  for(bin_var in bin_variables) {
    for(cont_var in cont_variables) {
      invisible(capture.output( aux <- psych::biserial(x = data_z[[cont_var]]
                                                , y = data[[bin_var]])))
      
      Sigma[cont_var, bin_var] <- Sigma[bin_var, cont_var] <- aux
      
    }
    for(cat_variables in categ_variables) {
      invisible(capture.output( aux <- psych::polychoric(
        x = data[, c(cat_variables, bin_var)],global = FALSE)$rho[1, 2]))
      
      Sigma[cat_variables, bin_var] <- Sigma[bin_var, cat_variables] <- aux
      
    }
    for(second_bin_var in bin_variables) {

      if(second_bin_var != bin_var) {
        invisible(capture.output( aux <- psych::tetrachoric(
          x = data[,c(bin_var,second_bin_var)])$rho[1, 2]))
        
        Sigma[second_bin_var, bin_var] <- Sigma[bin_var, second_bin_var] <- aux
      }
    }
  }
  #Categorical correlations
  for(cat_var in categ_variables){
    for(cont_var in cont_variables) {
      invisible(capture.output( aux <- psych::polyserial(x = data_z[[cont_var]],
                                                  y = as.matrix(data[[cat_var]]))))
      
      Sigma[cont_var, cat_var] <- Sigma[cat_var, cont_var] <- aux
      
    }
    for(second_cat_var in categ_variables) {
      if(second_cat_var != cat_var) {
        
        invisible(capture.output( aux <- psych::polychoric(
          x = data[, c(cat_var,second_cat_var)])$rho[1, 2]))
        
        Sigma[second_cat_var, cat_var] <- Sigma[cat_var, second_cat_var] <- aux
      }
    }
  }


return(Sigma)
}