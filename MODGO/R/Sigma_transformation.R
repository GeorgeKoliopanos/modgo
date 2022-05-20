#' Correlation transformation
#' 
#' Transforming covariance matrix using biserial, tetrachoric and
#' polyserial functions 
#' 
#' @param x a vector 
#' @param ties_method Method on how to deal with equal values 
#' during rank transformation.Acceptable input:"max","average","min".
#' @return A numeric vector.
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