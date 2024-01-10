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
#' @param n_samples Number of rows of each simulated data set. Default is
#' the number of rows of \code{data}.
#' @param gener_var_model A matrix that contains two columns named "Variable" and
#' "Model". This matrix can be used only if a gener_var_model argument is
#' provided. It specifies what model should be used for each Variable.
#' Model values should be "RMFMKL", "RPRS", "STAR" or a combination of them,
#' e.g. "RMFMKL-RPRS" or "STAR-STAR", in case the use wants a bimodal simulation.
#' The user can select Generalised Poisson model for poisson variabes,
#' but this model cannot be included in bimodal simulation
#' @param gener_var_lmbds A matrix that contains lmbds values for each of the
#' variables of the data set to be used for either Generalized Lambda Distribution
#' Generalized Poisson Distribution or setting up thresholds 
#' @return A correlation matrix.
#' @author Francisco M. Ojeda, George Koliopanos
#' @keywords Normal rank transformation



Inverse_transformation_variables <- function (data, df_sim,
                                              variables, bin_variables,
                                              categ_variables, n_samples,
                                              gener_var, gener_var_lmbds) {
  
  for (j in 1:length(variables)) {
    variable <- colnames(df_sim)[j]
    if (!(gener_var)) {
      # Default modgo rank transform
      df_sim[[j]] <- rbi_normal_transform_inv(df_sim[[j]], data[[j]])
      
    } else if (gener_var) {
      # Generalized Lambda simulation = TRUE
      df_sim[[j]] <- general_transform_inv(x = df_sim[[j]],
                                           data = data,
                                           n_samples = n_samples,
                                           lmbds = gener_var_lmbds[, variable])
      # Round categorical variables depending
      if (variable %in% categ_variables) {
        if (!(variable %in% count_variables)) {
          
          df_sim[[j]] <- round(df_sim[[j]])
          df_sim[[j]][df_sim[[j]] > max(data[[j]])] <-
            max(data[[j]])
          df_sim[[j]][df_sim[[j]] < min(data[[j]])] <-
            min(data[[j]])
        } else{
          if (gener_var_lmbds[1, variable] >= 10) {
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