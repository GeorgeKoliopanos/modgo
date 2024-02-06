#' Inverse transform variables
#'
#' This function is used internally by \code{\link[modgo]{modgo}}. It transforms
#' all variables to their original scale.
#'
#' @param data a data frame with original variables.
#' @param df_sim data frame with transformed variables.
#' @param variables variables a character vector indicating which
#' columns of \code{data} should be used.
#' @param bin_variables a character vector listing the binary variables.
#' @param  categ_variables a character vector listing the ordinal categorical
#' variables.
#' @param count_variables a character vector listing the count as a sub
#'  sub category of categorical variables. Count variables should be part
#'  of categorical variables vector. Count variables are treated differently
#'  when using gldex to simulate them.
#' @param n_samples Number of rows of each simulated data set. Default is
#' the number of rows of \code{data}.
#' @param generalized_mode A logical value indicating if 
#' generalized lambda/poisson distributions or set up thresholds will be used to
#' generate the simulated values
#' @param generalized_mode_lmbds A matrix that contains lambdas values for each 
#' of the variables of the data set to be used for either Generalized Lambda
#' Distribution Generalized Poisson Distribution or setting up thresholds 
#' @return A correlation matrix.
#' @author Francisco M. Ojeda, George Koliopanos
#' @keywords Normal rank transformation



Inverse_transformation_variables <- function (data,
                                              df_sim,
                                              variables,
                                              bin_variables,
                                              categ_variables,
                                              count_variables,
                                              n_samples,
                                              generalized_mode,
                                              generalized_mode_lmbds) {
  
  for (j in 1:length(variables)) {
    variable <- colnames(df_sim)[j]
    if (!(generalized_mode)) {
      # Default modgo rank transform
      df_sim[[j]] <- rbi_normal_transform_inv(df_sim[[j]], data[[j]])
      
    } else if (generalized_mode) {
      # Generalized Lambda simulation = TRUE
      df_sim[[j]] <- general_transform_inv(x = df_sim[[j]],
                                           data = data,
                                           n_samples = n_samples,
                                           lmbds = generalized_mode_lmbds[, variable])
      # Round categorical variables depending
      if (variable %in% categ_variables) {
        if (!(variable %in% count_variables)) {
          
          df_sim[[j]] <- round(df_sim[[j]])
          df_sim[[j]][df_sim[[j]] > max(data[[j]])] <-
            max(data[[j]])
          df_sim[[j]][df_sim[[j]] < min(data[[j]])] <-
            min(data[[j]])
        } else{
          if (generalized_mode_lmbds[1, variable] >= 10) {
            df_sim[[j]] <- round(df_sim[[j]])
            df_sim[[j]][df_sim[[j]] < 0] <- 0
          } else{
            df_sim[[j]] <- floor(df_sim[[j]])
            df_sim[[j]][df_sim[[j]] < 0] <- 0
          }
        }
        
      }
    }
  }
 
}