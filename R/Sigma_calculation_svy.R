#' Calculate Sigma with the help of wCorr::weightedCorr function
#'
#' This function is used internally by \code{\link[modgo]{modgo_svy}}. It 
#' conducts the computation of the correlation matrix of the 
#' transformed variables, which are assumed to follow a multivariate 
#' normal distribution.
#'
#' @param design survey design object from survey package containing the 
#' data to be used.
#' @param variables variables a character vector indicating which
#' variables of \code{design} should be used.
#' @param bin_variables a character vector listing the binary variables.
#' @param  categ_variables a character vector listing the ordinal categorical
#' variables.
#' @return A correlation matrix.
#' @author Francisco M. Ojeda, George Koliopanos
#' @keywords Normal rank transformation



Sigma_calculation_svy <- function (design,
                                   variables,
                                   bin_variables,
                                   categ_variables) {
  
  if (!requireNamespace("wCorr", quietly = TRUE)) {
    stop(
      "Package \"wCorr\" must be installed to use this function.",
      call. = FALSE
    )
  }   
  
  
  # Find the continuous variables
  continuous_var <-
    setdiff(variables, c(bin_variables, categ_variables))  
  
  # Rank transform original data
  
  design_rbi <- design
  for(i in variables) {
    design_rbi <- rbi_normal_transform_svy(
      x = i, design = design_rbi, x_transf = i
    )
  }
  

  # Calculate correlation matrix
  Sigma <- cov(design[["variables"]][, variables])
  diag(Sigma) <- 1
  
  # Check if categorical variable has
  # more than 7 values. If so it will be treated as continuous
  categ_var <- NULL
  for (ct_var in categ_variables) {
    if (length(unique(design[["variables"]][, ct_var])) > 7) {
      message(
        paste0(
          ct_var,
          " will be treated as continuous due to having more than 7 unique",
          "values."
        )
      )
    } else {
      categ_var <- c(categ_var, ct_var)
    }
  }
  
  # Calculate correlation matrix using wCorr::weightedCorr
  Sigma <- Sigma_transformation_svy(
    design = design,
    design_z = design_rbi,
    Sigma,
    variables = variables,
    bin_variables = bin_variables,
    categ_variables = categ_var
  )
  

  return(Sigma)
}