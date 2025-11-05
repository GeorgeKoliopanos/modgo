#' Correlation of transformed variables
#'
#' This function is used internally by \code{\link[modgo]{modgo_svy}}. It 
#' finishes the computation of the correlation matrix of the transformed 
#' variables, which
#' are assumed to follow a multivariate normal distribution. It computes the
#' correlations involving at least one categorical variable. For this purpose
#' the polyserial and polychoric correlations are used.
#'
#' @param design survey design object from survey package 
#' containing the original data.
#' @param design_z survey design object from survey package 
#' containing the transformed data.
#' @param Sigma A numeric square matrix.
#' @param variables variables a character vector indicating which
#' variables of \code{design} should be used.
#' @param bin_variables a character vector listing the binary variables.
#' @param  categ_variables a character vector listing the ordinal categorical
#' variables.
#' @return A correlation matrix.
#' @author Francisco M. Ojeda, George Koliopanos
#' @keywords Normal rank transformation
#'


Sigma_transformation_svy <- function(design,
                                     design_z,
                                     Sigma,
                                     variables,
                                     bin_variables = c(),
                                     categ_variables = c()) {
  
  data <- design[["variables"]]
  data_z <- design_z[["variables"]]
  w <- weights(design, type = "sampling")
  
  cont_variables <-
    setdiff(variables, c(bin_variables, categ_variables))
  # Binary and categoric correlations
  
  for (bin_cat_var in c(bin_variables, categ_variables)) {
    for (cont_var in cont_variables) {
      
      aux <- wCorr::weightedCorr(x = data_z[[cont_var]], 
                                 y = data[[bin_cat_var]],
                                 method = "Polyserial",
                                 weights = w)    
      Sigma[cont_var, bin_cat_var] <- Sigma[bin_cat_var, cont_var] <- aux
      
    }

    for (second_bin_cat_var in c(bin_variables, categ_variables)) {
      if (second_bin_cat_var != bin_cat_var) {
        aux <- wCorr::weightedCorr(x = data[[second_bin_cat_var]], 
                                   y = data[[bin_cat_var]],
                                   method = "Polychoric",
                                   weights = w)            
        Sigma[second_bin_cat_var, bin_cat_var] <-
          Sigma[bin_cat_var, second_bin_cat_var] <- aux
      }
    }
  }
  

  return(Sigma)
}