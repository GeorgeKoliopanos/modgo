#' Correlation of transformed variables
#'
#' This function is used internally by \code{\link[modgo]{modgo}}. It finishes
#' the computation of the correlation matrix of the transformed variables, which
#' are assumed to follow a multivariate normal distribution. It computes the
#' correlations involving at least one categorical variable. For this purpose
#' the biserial, tetrachoric, polyserial and polychoric correlations are used.
#'
#' @param data a data frame with original variables.
#' @param data_z data frame with transformed variables.
#' @param Sigma A numeric square matrix.
#' @param variables variables a character vector indicating which
#' columns of \code{data} should be used.
#' @param bin_variables a character vector listing the binary variables.
#' @param  categ_variables a character vector listing the ordinal categorical
#' variables.
#' @return A correlation matrix.
#' @author Francisco M. Ojeda, George Koliopanos
#' @keywords Normal rank transformation
#'
#'@importFrom psych biserial
#'@importFrom psych polychoric
#'@importFrom psych tetrachoric
#'@importFrom psych polyserial
#'@importFrom utils capture.output

Sigma_transformation <- function(data,
                                 data_z,
                                 Sigma,
                                 variables,
                                 bin_variables = c(),
                                 categ_variables = c()) {
  cont_variables <-
    setdiff(variables, c(bin_variables, categ_variables))
  #Binary correlations
  for (bin_var in bin_variables) {
    for (cont_var in cont_variables) {
      invisible(capture.output(aux <-
                                 psych::biserial(x = data_z[[cont_var]]
                                                 , y = data[[bin_var]])))
      
      Sigma[cont_var, bin_var] <- Sigma[bin_var, cont_var] <- aux
      
    }
    for (cat_variables in categ_variables) {
      invisible(capture.output(aux <- psych::polychoric(x = data[, c(cat_variables, bin_var)], global = FALSE)$rho[1, 2]))
      
      Sigma[cat_variables, bin_var] <-
        Sigma[bin_var, cat_variables] <- aux
      
    }
    for (second_bin_var in bin_variables) {
      if (second_bin_var != bin_var) {
        invisible(capture.output(aux <- psych::tetrachoric(x = data[, c(bin_var, second_bin_var)])$rho[1, 2]))
        
        Sigma[second_bin_var, bin_var] <-
          Sigma[bin_var, second_bin_var] <- aux
      }
    }
  }
  #Categorical correlations
  for (cat_var in categ_variables) {
    for (cont_var in cont_variables) {
      invisible(capture.output(aux <-
                                 psych::polyserial(
                                   x = data_z[[cont_var]],
                                   y = as.matrix(data[[cat_var]])
                                 )))
      
      Sigma[cont_var, cat_var] <- Sigma[cat_var, cont_var] <- aux
      
    }
    for (second_cat_var in categ_variables) {
      if (second_cat_var != cat_var) {
        invisible(capture.output(aux <- psych::polychoric(x = data[, c(cat_var, second_cat_var)])$rho[1, 2]))
        
        Sigma[second_cat_var, cat_var] <-
          Sigma[cat_var, second_cat_var] <- aux
      }
    }
  }
  
  
  return(Sigma)
}