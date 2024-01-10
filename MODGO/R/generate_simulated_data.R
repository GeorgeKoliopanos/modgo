#' Generate new data set by using previous correlation matrix
#'
#' This function is used internally by \code{\link[modgo]{modgo}}. It conducts
#' the computation of the correlation matrix of the transformed variables, which
#' are assumed to follow a multivariate normal distribution.
#'
#' @param data a data frame with original variables.
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
#' @param gener_var A logical value indicating if generalized lambda/poisson
#'  distributions or set up thresholds will be used to generate the simulated values
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
#' @param multi_sugg_prop A named vector that provides a  proportion of
#'  value=1 for specific binary variables(=name of the vector) that will be
#'  the close to the proportion of this value in the simulated data sets.
#' @param pertr_vec A named vector.Vector's names are the continuous variables
#' that the user want to perturb. Variance of simulated data set mimic original
#' data's variance.
#' @param var_infl A named vector.Vector's names are the continuous variables
#' that the user want to perturb and increase their variance
#' @param infl_cov_stable Logical value. If TRUE,perturbation is applied to
#' original data set and simulations values mimic the perturbed original data
#' set.Covariance matrix used for simulation = original data's correlations.
#' If FALSE, perturbation is applied to the simulated data sets.
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
                           gener_var,
                           gener_var_lmbds,
                           multi_sugg_prop,
                           pertr_vec,
                           var_infl,
                           infl_cov_stable
                          ){
  for (j in 1:length(variables)) {
    variable <- variables[[j]]
    # Default modgo rank inverse transformation
    if(gener_var == FALSE){
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
      } else if (gener_var == TRUE) {
        # Generalised transformation using Generalised Lambdas
        df_sim[[j]] <- general_transform_inv(x = df_sim[[j]],
                                             data = data[[j]],
                                             n_samples = n_samples,
                                             lmbds = gener_var_lmbds[, variable])
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
            if (gener_var_lmbds[1, variable] >= 10) {
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