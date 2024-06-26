#' Generate new dataset by using previous correlation matrix
#'
#' This function is used internally by \code{\link[modgo]{modgo}}. It conducts
#' the computation of the correlation matrix of the transformed variables, which
#' are assumed to follow a multivariate normal distribution.
#'
#' @param data a data frame with original variables.
#' @param df_sim a data frame with simulated values.
#' @param variables variables a character vector indicating which
#' columns of \code{data} should be used.
#' @param bin_variables a character vector listing the binary variables.
#' @param  categ_variables a character vector listing the ordinal categorical
#' variables.
#' @param count_variables a character vector listing the count as a sub
#'  sub category of categorical variables. Count variables should be part
#'  of categorical variables vector. Count variables are treated differently
#'  when using gldex to simulate them.
#' @param n_samples Number of rows of each simulated dataset. Default is
#' the number of rows of \code{data}.
#' @param generalized_mode A logical value indicating if generalized lambda/Poisson
#'  distributions or set up thresholds will be used to generate the simulated values
#' @param generalized_mode_lmbds A matrix that contains lmbds values for each of the
#' variables of the dataset to be used for either Generalized Lambda Distribution
#' Generalized Poisson Distribution or setting up thresholds
#' @param multi_sugg_prop A named vector that provides a  proportion of
#'  value=1 for specific binary variables(=name of the vector) that will be
#'  the close to the proportion of this value in the simulated datasets.
#' @param pertr_vec A named vector.Vector's names are the continuous variables
#' that the user want to perturb. Variance of simulated dataset mimic original
#' data's variance.
#' @param var_infl A named vector.Vector's names are the continuous variables
#' that the user want to perturb and increase their variance
#' @param infl_cov_stable Logical value. If TRUE,perturbation is applied to
#' original dataset and simulations values mimic the perturbed original data
#' set.Covariance matrix used for simulation = original data's correlations.
#' If FALSE, perturbation is applied to the simulated datasets.
#' @return Simulation Data Frame
#' @author Francisco M. Ojeda, George Koliopanos
#' @keywords Normal rank transformation

generate_simulated_data <- function(data,
                           df_sim,
                           variables,
                           bin_variables,
                           categ_variables,
                           count_variables,
                           n_samples,
                           generalized_mode,
                           generalized_mode_lmbds,
                           multi_sugg_prop,
                           pertr_vec,
                           var_infl,
                           infl_cov_stable
                          ){
  for (j in 1:length(variables)) {
    variable <- variables[[j]]
    # Default modgo rank inverse transformation
    if(generalized_mode == FALSE){
      # Multi suggestive proportion rank inverse transormation 
      if (colnames(df_sim)[[j]] %in% names(multi_sugg_prop)) {
        df_sim[[j]] <- rbi_normal_transform_inv(df_sim[[j]],
                                                rbinom(n = dim(df_sim)[1],
                                                       1,
                                                       prob = multi_sugg_prop[variable]))
      } else {
        # Default rank inverse transformation
          df_sim[[j]] <- rbi_normal_transform_inv(df_sim[[j]],
                                                  data[[j]])
        } 
      } else if (generalized_mode == TRUE) {
        # Generalised transformation using Generalised Lambdas
        df_sim[[j]] <- general_transform_inv(x = df_sim[[j]],
                                             data = data[[j]],
                                             n_samples = n_samples,
                                             lmbds = generalized_mode_lmbds[, variable])
      }
    # Round categorical simulated data
    if (variable %in% categ_variables) {
          if (!(variable %in% count_variables)) {
            df_sim[[j]] <- round(df_sim[[j]])
            if (!is.null(data)){
              df_sim[[j]][df_sim[[j]] > max(data[[j]])] <-
                max(data[[j]])
              df_sim[[j]][df_sim[[j]] < min(data[[j]])] <-
                min(data[[j]])
            }
          } else {
            # Round/Floor count variables depending on theta value
            if (generalized_mode_lmbds[1, variable] >= 10) {
              df_sim[[j]] <- round(df_sim[[j]])
              df_sim[[j]][df_sim[[j]] < 0] <- 0
            } else{
              df_sim[[j]] <- floor(df_sim[[j]])
              df_sim[[j]][df_sim[[j]] < 0] <- 0
            }
          }
        }
      
    # Perturbation of simulated values
    if (variables[[j]] %in% names(pertr_vec)) {
      p <- pertr_vec[which(names(pertr_vec) == variables[[j]])]
      
      df_sim[[j]] <- (df_sim[[j]] * sqrt(1 - p)) +
        rnorm(length(df_sim[[j]]),
              mean = 0,
              sd = sd(df_sim[[j]]) * sqrt(p))
    }
    # Inflation of simulated data
    if (variables[[j]] %in% names(var_infl) &&
        infl_cov_stable == FALSE) {
      p <- var_infl[which(names(var_infl) == variables[[j]])]
      
      df_sim[[j]] <- df_sim[[j]] +
        rnorm(length(df_sim[[j]]),
              mean = 0,
              sd = sd(sqrt(p) * df_sim[[j]]))
    }
  }
  return(df_sim)
}