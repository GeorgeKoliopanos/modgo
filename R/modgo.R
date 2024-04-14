#' MOck Data GeneratiOn
#'
#' Creates synthetic dataset based on real data by means of the
#' rank-based inverse normal transformation. Data with perturbed
#' characteristics can be generated.
#'
#' Simulated data is generated based on available data. The algorithm 
#' used is based on the ranked-based inverse normal transformation 
#' (Koliopanos et al. (2023)) and attempts to mimic the characteristics of the 
#' original data. 
#' 
#' All variables passed to \code{modgo} should be of class 
#' double or integer. This includes the variables passed to the parameter
#' \code{categ_variables}. The character vector \code{variables}, indicating
#' the variables in \code{data} to be used in the simulation, should 
#' contain at least two variables. The variables in \code{variables} not present
#' in \code{bin_variables} nor \code{categ_variables} will be treated as 
#' continuous variables.
#'
#' @param data A data frame containing the data whose characteristics are to be
#' mimicked during the data simulation.
#' @param sigma A covariance matrix of NxN (N= number of variables)
#' provided by the user to bypass the covariance matrix calculations
#' @param ties_method Method used to deal with ties
#' during rank transformation. Allowed input: "max","average" or "min". This
#' parameter is passed by \code{\link[modgo]{rbi_normal_transform}} to the
#' parameter \code{ties.method} of \code{\link[base]{rank}}.
#' @param variables A character vector indicating the columns in \code{data} 
#' to be used. Default: \code{colnames(data)}.
#' @param bin_variables  A character vector listing those entries in 
#' \code{variables} to be treated as binary variables. 
#' @param categ_variables A character vector listing those entries in 
#' \code{variables} to be treated as ordinal categorical variables, with 
#' more than two categories. See Details.
#' @param count_variables A character vector listing those entries 
#' \code{categ_variables}  to be treated as count variables. Relevant only when
#' \code{generalized_mode = TRUE}.
#' @param nrep Number of simulated datasets to be generated.
#' @param noise_mu Logical. Should noise be added to the  
#' mean vector of the multivariate normal distribution used to draw the 
#' simulated values? Default: FALSE.
#' @param pertr_vec A named vector. Vector's names are the continuous variables
#' that the user want to perturb. Variance of simulated dataset mimic original
#' data's variance.
#' @param var_infl A named vector. Vector's names are the continuous variables
#' that the user want to perturb and increase their variance
#' @param infl_cov_stable Logical value. If TRUE,perturbation is applied to
#' original dataset and simulations values mimic the perturbed original 
#' dataset. Covariance matrix used for simulation = original data's correlations.
#' If FALSE, perturbation is applied to the simulated datasets.
#' @param n_samples Number of rows of each simulated dataset. Default is
#' the number of rows of \code{data}.
#' @param change_cov Change the covariance of a specific pair of variables.
#' @param change_amount the amount of change in  the covariance
#'  of a specific pair of variables.
#' @param seed A numeric value specifying the random seed. If \code{seed = NA},
#' no random seed is set.
#' @param thresh_var A data frame that contains the thresholds(left and right)
#' of specified variables
#' (1st column: variable names, 2nd column: Left thresholds,
#' 3rd column: Right thresholds)
#' @param thresh_force A logical value indicating if you want to force threshold
#' in case the proportion of samples that can surpass the threshold are less
#' than 10\%
#' @param var_prop A named vector that provides a  proportion of
#'  value=1 for a specific binary variable (=name of the vector) that will be
#'  the proportion of this value in the simulated datasets.[this may increase
#'  execution time drastically]
#' @param multi_sugg_prop A named vector that provides a  proportion of
#'  value=1 for specific binary variables (=name of the vector) that will be
#'  the close to the proportion of this value in the simulated datasets.
#' @param tol A numeric value that set up
#'  tolerance(relative to largest variance) for numerical lack of
#'  positive-definiteness in Sigma
#' @param stop_sim A logical value indicating if the analysis should
#' stop before simulation and produce only the correlation matrix
#' @param generalized_mode A logical value indicating if generalized lambda/Poisson
#'  distributions or set up thresholds will be used to generate the simulated values
#' @param generalized_mode_model A matrix that contains two columns named 
#' "Variable" and "Model". This matrix can be used only if a 
#' generalized_mode_model argument is provided. It specifies what model should 
#' be used for each Variable. Model values should be "rmfmkl", "rprs", "star" 
#' or a combination of them, e.g. "rmfmkl-rprs" or "star-star", in case the 
#' user wants a bimodal simulation.
#' The user can select Generalised Poisson model for Poisson variables,
#' but this model cannot be included in bimodal simulation
#' @param generalized_mode_lmbds A matrix that contains lambdas values for each 
#' of the variables of the dataset to be used for either Generalized Lambda 
#' Distribution Generalized Poisson Distribution or setting up thresholds
#' @param new_mean_sd A matrix that contains two columns named
#' "Mean" and "SD" that the user specifies desired Means and Standard Deviations
#' in the simulated datasets for specific continues variables. The variables
#' must be declared as ROWNAMES in the matrix.
#' @return A list with the following components:
#' \item{simulated_data}{A list of data frames containing the simulated data.}
#' \item{original_data}{A data frame with the input data.}
#' \item{correlations}{A list of correlation matrices. The ith element is the
#' correlation matrix for the ith simulated dataset. The \code{(repn + 1)}the
#' (last) element of the list is the average of the correlation matrices.}
#' \item{bin_variables}{A character vector listing the binary variables}
#' \item{categ_variables }{A character vector listing the ordinal
#' categorical variables}
#' \item{covariance_matrix}{Covariance matrix used when generating observations
#' from a multivariate normal distribution.}
#' \item{seed}{Random seed used.}
#' \item{samples_produced}{Number of rows of each simulated dataset.}
#' \item{sim_dataset_number}{Number of simulated datasets produced.}
#' @author Francisco M. Ojeda, George Koliopanos
#' @keywords mock data generation
#' @references Koliopanos, G., Ojeda, F. and Ziegler A. (2023).
#' A simple-to-use R package for mimicking study data by simulations.
#' \emph{Methods Inf Med}, 62(03/04), 119-129.
#' @examples
#' data("Cleveland",package="modgo")
#' test_modgo <- modgo(data = Cleveland,
#'      bin_variables = c("CAD","HighFastBloodSugar","Sex","ExInducedAngina"),
#'      categ_variables =c("Chestpaintype"))
#' @export
#' @importFrom Matrix nearPD
#' @importFrom MASS mvrnorm
#' @import stats

modgo <-
  function(data,
           ties_method =  "max",
           variables = colnames(data),
           bin_variables = NULL,
           categ_variables = NULL,
           count_variables = NULL,
           n_samples = nrow(data),
           sigma = NULL,
           nrep = 100,
           noise_mu = FALSE,
           pertr_vec = NULL,
           change_cov = NULL,
           change_amount = 0,
           seed = 1,
           thresh_var = NULL,
           thresh_force = FALSE,
           var_prop = NULL,
           var_infl = NULL,
           infl_cov_stable = FALSE,
           tol = 1e-06,
           stop_sim = FALSE,
           new_mean_sd = NULL,
           multi_sugg_prop = NULL,
           generalized_mode = FALSE,
           generalized_mode_model = NULL,
           generalized_mode_lmbds = NULL) {
    if (!is.na(seed)) {
      # Setting Seed
      set.seed(seed)
      
    }

    # Check Arguments
    .args <- as.list(match.call()[-1])
    do.call(checkArguments, .args, envir = parent.frame(n = 2))

    # Find the continuous variables
    continuous_var <-
      setdiff(variables, c(bin_variables, categ_variables))
    
    data <- data[, variables]        
    
    # Noise in multivariate distributions centers
    if (noise_mu == TRUE) {
      ns <- rnorm(length(variables), mean = 0, sd = 1)
      
      
    } else {
      ns <- matrix(0, nrow = 1, ncol = length(variables))
      
    }
    
    OriginalData <- data
    
    # Normal run in case user does not provide a sigma
    if (length(sigma) == 0) {
      # Calculate correlation matrix with the use of polychoric/polyserial
      Sigma <- Sigma_calculation(data = data,
                                 variables = variables,
                                 bin_variables = bin_variables,
                                 categ_variables = categ_variables,
                                 ties_method = ties_method)
    } else {
      # Set sigma if provided by user
      Sigma <- sigma
    }
    # Save sigma before near positive definite matrix and set of specified pairs
    
    pre_sim_sigma <- Sigma
    
    # Change specified pairs of covariance matrix by specified amount
    if (change_amount != 0 && length(change_cov) == 2) {
      Sigma[change_cov[1], change_cov[2]] <-
        Sigma[change_cov[2], change_cov[1]] <-
        Sigma[change_cov[2], change_cov[1]] + change_amount
    }
    
    
    
    # Check Sigma if it is positive definite
    if (!all(eigen(Sigma)$value > 0)) {
      message("Covariance matrix is not positive definite.")
      message("It will be replaced with nearest positive definite matrix")
      Sigma <- Matrix::nearPD(Sigma, corr = TRUE)$mat
    }
    
    # Check Sigma if it contains NA's
    if (any(is.na(Sigma))) {
      stop("Correlation of data contains NA's")
    }
    
    
    ## Inflation analysis - stable covariance matrix
    if (length(var_infl) > 0 && infl_cov_stable == TRUE) {
      #Inverse transformation of each variable
      for (j in 1:length(variables)) {
        if (colnames(data)[[j]] %in% names(var_infl)) {
          p <- var_infl[which(names(var_infl) == colnames(data)[[j]])]
          data[[j]] <- data[[j]] +
            rnorm(n_samples,
                  mean = 0,
                  sd = sd(sqrt(p) * data[[j]]))
          
        }
      }
    }
    
    # Prepare generalized_mode four moments
    if (is.null(generalized_mode_lmbds) && generalized_mode == TRUE) {
      generalized_mode_lmbds <- generalizedMatrix(data,
                                           variables,
                                           bin_variables,
                                           generalized_mode_model,
                                           multi_sugg_prop)
    }
    
    ## Thresholds
    if (length(thresh_var[, 1]) > 0) {
      mt_sim <- MASS::mvrnorm(
        n = n_samples,
        mu = rep(0, length(variables)) + ns,
        Sigma = Sigma,
        tol = tol
      )
      df_sim <- data.frame(mt_sim)
      names(df_sim) <- variables
      #Inverse transformation of each variable
      df_sim <- generate_simulated_data(data = data,
                              df_sim = df_sim,
                              variables = variables,
                              bin_variables = bin_variables,
                              categ_variables = categ_variables,
                              count_variables = count_variables,
                              n_samples = n_samples,
                              generalized_mode = generalized_mode,
                              generalized_mode_lmbds = generalized_mode_lmbds,
                              multi_sugg_prop = NULL,
                              pertr_vec = NULL,
                              var_infl = NULL,
                              infl_cov_stable = NULL)
      
      for (i in c(1:length(thresh_var[, 1]))) {
        if (is.na(thresh_var[i, 2]))
        {
          low_thresh <-
            min(df_sim[thresh_var[i, 1]]) - 1
        } else{
          low_thresh <- thresh_var[i, 2]
        }
        if (is.na(thresh_var[i, 3]))
        {
          up_thresh <-
            max(df_sim[thresh_var[i, 1]]) + 1
        } else{
          up_thresh <- thresh_var[i, 3]
        }
        
        df_sim <-
          df_sim[which(df_sim[thresh_var[i, 1]] < up_thresh &
                         df_sim[thresh_var[i, 1]] > low_thresh),]
      }
      
      if (length(df_sim[, 1]) / n_samples < 0.1 &
          thresh_force == FALSE) {
        stop("The propotion of simulated samples passing
           the threshold is less than 10% ")
      } else if (length(df_sim[, 1]) == 0) {
        stop("The propotion of simulated samples passing
           the threshold is 0% ")
      } else{
        thresh_multi <- 1.1 * n_samples / length(df_sim[, 1])
      }
      
    } else{
      thresh_multi = 1
    }
    
    # Binary variables proportions
    if (length(var_prop) == 1) {
      mt_sim <- MASS::mvrnorm(
        n = n_samples,
        mu = rep(0, length(variables)) + ns,
        Sigma = Sigma,
        tol = tol
      )
      df_sim <- data.frame(mt_sim)
      names(df_sim) <- variables
      #Inverse transformation of each variable
      df_sim <- generate_simulated_data(data = data,
                              df_sim = df_sim,
                              variables = variables,
                              bin_variables = bin_variables,
                              categ_variables = categ_variables,
                              count_variables = count_variables,
                              n_samples = n_samples,
                              generalized_mode = generalized_mode,
                              generalized_mode_lmbds = generalized_mode_lmbds,
                              multi_sugg_prop = NULL,
                              pertr_vec = NULL,
                              var_infl = NULL,
                              infl_cov_stable = NULL)
      
      counts_1 <- length(which(df_sim[, names(var_prop)] == 1))
      counts_0 <- length(which(df_sim[, names(var_prop)] == 0))
      req_1 <- var_prop * n_samples
      req_0 <- (1 - var_prop) * n_samples
      
      if (req_1 >= counts_1) {
        thresh_multi <-  (1.1 * req_1) / counts_1
      } else if (req_1 < counts_1) {
        thresh_multi <-  (1.1 * req_0) / counts_0
      }
    }
    
    # Starting loop for many repetions
    
    Correlations <- vector(mode = "list", length = nrep + 1)
    mean_corr <- matrix(0, nrow = length(variables), ncol = length(variables))
    
    SimulatedData <- vector(mode = "list", length = nrep)
    i <- 1
    counter <- 0
    if (stop_sim == FALSE) {
      # Loop for creating new datasets and obtaining their mean correlations
      while (i < nrep + 1) {
        mt_sim <- MASS::mvrnorm(
          n = ceiling(n_samples * thresh_multi),
          mu = rep(0, length(variables)) + ns,
          Sigma = Sigma,
          tol = tol
        )
        df_sim <- data.frame(mt_sim)
        names(df_sim) <- variables
        #Inverse transformation of each variable
        df_sim <- generate_simulated_data(data = data,
                                df_sim = df_sim,
                                variables = variables,
                                bin_variables = bin_variables,
                                categ_variables = categ_variables,
                                count_variables = count_variables,
                                n_samples = n_samples,
                                generalized_mode = generalized_mode,
                                generalized_mode_lmbds = generalized_mode_lmbds,
                                multi_sugg_prop = multi_sugg_prop,
                                pertr_vec = pertr_vec,
                                var_infl = var_infl,
                                infl_cov_stable = infl_cov_stable)
        #Threshold process
        if (length(thresh_var[, 1]) > 0) {
          for (j in c(1:length(thresh_var[, 1]))) {
            if (is.na(thresh_var[j, 2]))
            {
              low_thresh <-
                min(df_sim[thresh_var[j, 1]]) - 1
            } else{
              low_thresh <- thresh_var[j, 2]
            }
            if (is.na(thresh_var[j, 3]))
            {
              up_thresh <-
                max(df_sim[thresh_var[j, 1]]) + 1
            } else{
              up_thresh <- thresh_var[j, 3]
            }
            df_sim <-
              df_sim[which(df_sim[thresh_var[j, 1]] < up_thresh &
                             df_sim[thresh_var[j, 1]] > low_thresh),]
          }
          
          
        }
        #Proportion process
        
        if (length(var_prop) == 1) {
          rounded_length <- ceiling(n_samples * var_prop)
          df_sim_1 <-
            df_sim[which(df_sim[, names(var_prop)] == 1),]
          df_sim_0 <-
            df_sim[which(df_sim[, names(var_prop)] == 0),]
          df_sim <- rbind(df_sim_1[c(1:rounded_length),],
                          df_sim_0[c(1:(n_samples - rounded_length)),])
          
        }
        if (length(df_sim[, 1]) < n_samples | is.null(df_sim) |
            any(apply(df_sim, 2, function(x)
              any(is.na(x)))) |
            any(is.na(suppressWarnings(cor(df_sim))))) {
          i <- i - 1
          
          counter <- counter + 1
        } else {
          if (!is.null(data)){
          df_sim <- df_sim[c(1:n_samples),]
          for (j in rownames(new_mean_sd)) {
            df_sim[[j]] <-
              ((df_sim[[j]] - mean(OriginalData[[j]])) /  sd(OriginalData[[j]])) * new_mean_sd[j, "SD"] +
              new_mean_sd[j, "Mean"]
          }
          }
          #Correlation calculation
          Correlations[[i]] <- cor(df_sim)
          SimulatedData[[i]] <- df_sim
          #Mean correlation calculation
          mean_corr <- mean_corr + (Correlations[[i]] / nrep)
        }
        
        if (counter > 100 * nrep) {
          stop("Could not create simulated datasets with this specific thresholds")
        }
        
        i <- i + 1
      }
      
      
      Correlations[[nrep + 1]] <- mean_corr
      names(Correlations) <- c(paste0("rep", seq(1:nrep)), "Mean")
      
      results <- list(
        SimulatedData,
        OriginalData,
        Correlations,
        bin_variables,
        categ_variables,
        Sigma,
        seed,
        n_samples,
        nrep,
        pre_sim_sigma,
        generalized_mode_lmbds
      )
      names(results) <- c(
        "simulated_data",
        "original_data",
        "correlations",
        "bin_variables",
        "categ_variables",
        "covariance_matrix",
        "seed",
        "samples_produced",
        "sim_dataset_number",
        "presim_sigma",
        "generalized_mode_lmbds"
      )
    } else{
      results <- list(
        OriginalData,
        bin_variables,
        categ_variables,
        Sigma,
        seed,
        pre_sim_sigma,
        generalized_mode_lmbds
      )
      names(results) <-
        c(
          "original_data",
          "bin_variables",
          "categ_variables",
          "covariance_matrix",
          "seed",
          "presim_sigma",
          "generalized_mode_lmbds"
        )
      
    }
    
    return(results)
  }
