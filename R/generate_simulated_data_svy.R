#' Invert rank-based inverse normal transformation to obtain simulated data
#'
#' This function is used internally by \code{\link[modgo]{modgo_svy}}. It 
#' takes data simulated from a multivariate normal distribution and transforms 
#' it back to the scale of the original data.
#'
#' @param design survey design object from survey package containing the 
#' original data.
#' @param df_sim a data frame with simulated values from a multivariate normal 
#' distribution
#' @param variables variables a character vector indicating which
#' columns of \code{data} should be used.
#' @param bin_variables a character vector listing the binary variables.
#' @param  categ_variables a character vector listing the ordinal categorical
#' variables.
#' @param n_samples Number of rows of each simulated dataset. Default is
#' the number of rows of \code{data}.
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

generate_simulated_data_svy <- function(design,
                                        df_sim,
                                        variables,
                                        bin_variables,
                                        categ_variables,
                                        n_samples,
                                        multi_sugg_prop,
                                        pertr_vec,
                                        var_infl,
                                        infl_cov_stable
){
  for (j in 1:length(variables)) {
    variable <- variables[[j]]
    # Default modgo rank inverse transformation

    # Multi suggestive proportion rank inverse transormation 
    if (variable %in% names(multi_sugg_prop)) {
      # USING REGULAR rbi_normal_transform_inv !!!!
      df_sim[[variable]] <- rbi_normal_transform_inv(
        x = df_sim[[variable]],
        rbinom(n = nrow(df_sim), 1, prob = multi_sugg_prop[variable])
      )
    } else {
      # Default rank inverse transformation
      df_sim[[variable]] <- rbi_normal_transform_inv_svy(
        x = df_sim[[variable]],
        design = design, 
        x_original = variable)
    } 
    # Round categorical simulated data
    if (variable %in% categ_variables) {
      
      df_sim[[variable]] <- round(df_sim[[variable]])
      if (!is.null(design)){
        max_tmp <- max(design[["variables"]][[variable]])
        min_tmp <- min(design[["variables"]][[variable]])
        df_sim[[variable]][df_sim[[variable]] > max_tmp] <- max_tmp
        df_sim[[j]][df_sim[[variable]] < min_tmp] <- min_tmp
        rm(max_tmp, min_tmp)
      }
      
    }
    
    # Perturbation of simulated values
    if (variable %in% names(pertr_vec)) {
      p <- pertr_vec[which(names(pertr_vec) == variable)]
      
      df_sim[[variable]] <- (df_sim[[variable]] * sqrt(1 - p)) +
        rnorm(length(df_sim[[variable]]),
              mean = 0,
              sd = sd(df_sim[[variable]]) * sqrt(p))
    }
    # Inflation of simulated data
    if (variable %in% names(var_infl) &&
        infl_cov_stable == FALSE) {
      p <- var_infl[which(names(var_infl) == variable)]
      
      df_sim[[variable]] <- df_sim[[variable]] +
        rnorm(length(df_sim[[variable]]),
              mean = 0,
              sd = sd(sqrt(p) * df_sim[[variable]]))
    }
  }
  return(df_sim)
}