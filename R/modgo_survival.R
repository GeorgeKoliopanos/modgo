#' MOck Data GeneratiOn
#'
#' \code{modgo_survival} Create mock dataset from a real one by using
#' Generalized Lambdas Distributions and by seperating the data set in 2 based
#' in the event status.
#'
#' Simulated data is generated based on available data. The simulated data
#' mimics the characteristics of the original data. The algorithm used is
#' based on the ranked based inverse normal transformation (Koliopanos et
#' al. (2023)).
#'
#' @param data a data frame containing the data whose characteristics are to be
#' mimicked during the data simulation.
#' @param sigma a covariance matrix of NxN (N= number of variables)
#' provided by the user to bypass the covariance matrix calculations
#' @param ties_method Method on how to deal with equal values
#' during rank transformation. Acceptable input:"max","average","min". This
#' parameter is passed by \code{\link[modgo]{rbi_normal_transform}} to the
#' parameter \code{ties.method} of \code{\link[base]{rank}}.
#' @param variables a vector of which variables you want to transform.
#' Default:colnames(data)
#' @param event_variable  a character string listing the event variable.
#' @param time_variable  a character string listing the time variable.
#' @param bin_variables  a character vector listing the binary variables.
#' @param categ_variables a character vector listing the ordinal categorical
#' variables.
#' @param count_variables a character vector listing the count as a sub
#'  sub category of categorical variables. Count variables should be part
#'  of categorical variables vector. Count variables are treated differently
#'  when using gldex to simulate them.
#' @param nrep number of repetitions.
#' @param noise_mu Logical value if you want to apply noise to
#' multivariate mean. Default: FALSE
#' @param pertr_vec A named vector.Vector's names are the continuous variables
#' that the user want to perturb. Variance of simulated data set mimic original
#' data's variance.
#' @param var_infl A named vector.Vector's names are the continuous variables
#' that the user want to perturb and increase their variance
#' @param infl_cov_stable Logical value. If TRUE,perturbation is applied to
#' original data set and simulations values mimic the perturbed original data
#' set.Covariance matrix used for simulation = original data's correlations.
#' If FALSE, perturbation is applied to the simulated data sets.
#' @param n_samples Number of rows of each simulated data set. Default is
#' the number of rows of \code{data}.
#' @param change_cov change the covariance of a specific pair of variables.
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
#'  value=1 for a specific binary variable(=name of the vector) that will be
#'  the proportion of this value in the simulated data sets.[this may increase
#'  execution time drastically]
#' @param multi_sugg_prop A named vector that provides a  proportion of
#'  value=1 for specific binary variables(=name of the vector) that will be
#'  the close to the proportion of this value in the simulated data sets.
#' @param tol A numeric value that set up
#'  tolerance(relative to largest variance) for numerical lack of
#'  positive-definiteness in Sigma
#' @param stop_sim A logical value indicating if the analysis should
#' stop before simulation and produce only the correlation matrix
#' @param generalized_mode A logical value indicating if generalized lambda/poisson
#'  distributions or set up thresholds will be used to generate the simulated values
#' @param generalized_mode_model A matrix that contains two columns named "Variable" and
#' "Model". This matrix can be used only if a generalized_mode_model argument is
#' provided. It specifies what model should be used for each Variable.
#' Model values should be "rmfmkl", "rprs", "star" or a combination of them,
#' e.g. "rmfmkl-rprs" or "star-star", in case the use wants a bimodal simulation.
#' The user can select Generalised Poisson model for poisson variables,
#' but this model cannot be included in bimodal simulation
#' @param generalized_mode_model_no_event A matrix that contains two columns named "Variable" and
#' "Model" and it is to be used for the non-event data set(event = 0). This matrix can be used only if a generalized_mode_model argument is
#' provided. It specifies what model should be used for each Variable.
#' Model values should be "rmfmkl", "rprs", "star" or a combination of them,
#' e.g. "rmfmkl-rprs" or "star-star", in case the use wants a bimodal simulation.
#' The user can select Generalised Poisson model for poisson variables,
#' but this model cannot be included in bimodal simulation
#' @param generalized_mode_model_event A matrix that contains two columns named "Variable" and
#' "Model" and it is to be used for the event data set(event = 1). This matrix can be used only if a generalized_mode_model argument is
#' provided. It specifies what model should be used for each Variable.
#' Model values should be "rmfmkl", "rprs", "star" or a combination of them,
#' e.g. "rmfmkl-rprs" or "star-star", in case the use wants a bimodal simulation.
#' The user can select Generalised Poisson model for poisson variables,
#' but this model cannot be included in bimodal simulation
#' @param generalized_mode_lmbds A matrix that contains lambdas values for each of the
#' variables of the data set to be used for either Generalized Lambda Distribution
#' Generalized Poisson Distribution or setting up thresholds
#' @param new_mean_sd A matrix that contains two columns named
#' "Mean" and "SD" that the user specifies desired Means and Standard Deviations
#' in the simulated data sets for specific continues variables. The variables
#' must be declared as ROWNAMES in the matrix
#' @param surv_method A numeric value that indicates which one of the 2 survival
#' methods will be used.
#' First method(surv_method = 1): Event and no event data sets are using 
#' different covariance matrices for the simulation.
#' Second method(surv_method = 2): Event and no event data sets
#' are using the same covariance matrix for the simulation 
#' @return A list with the following components:
#' \item{simulated_data}{A list of data frames containing the simulated data.}
#' \item{original_data}{A data frame with the input data.}
#' \item{correlations}{a list of correlation matrices. The ith element is the
#' correlation matrix for the ith simulated dataset. The \code{(repn + 1)}the
#' (last) element of the list is the average of the correlation matrices.}
#' \item{bin_variables}{character vector listing the binary variables}
#' \item{categ_variables }{a character vector listing the ordinal
#' categorical variables}
#' \item{covariance_matrix}{Covariance matrix used when generating observations
#' from a multivariate normal distribution.}
#' \item{seed}{Random seed used.}
#' \item{samples_produced}{Number of rows of each simulated dataset.}
#' \item{sim_dataset_number}{Number of simulated datasets produced.}
#' @author Francisco M. Ojeda, George Koliopanos
#' @keywords mock data generation
#' @export
#' @importFrom Matrix nearPD

modgo_survival <-
  function(data,
           event_variable = NULL,
           time_variable = NULL,
           surv_method = 1,
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
           generalized_mode = TRUE,
           generalized_mode_model = NULL,
           generalized_mode_model_event = "rprs",
           generalized_mode_model_no_event = "rprs",
           generalized_mode_lmbds = NULL) {
    
    .args <- as.list(match.call()[-1])
    # Check Arguments
    if (is.null(data)){
      stop("Data set is not provided")
    }
    
    # Event Variable checks
    if (length(event_variable) > 1) {
      stop("You can only specify one Event variable")
    }
    if (!(event_variable %in% variables)) {
      stop("Event variable is not part of the provided variables")
    }
    if (!(event_variable %in% bin_variables)) {
      stop("Event variable should be a binary variable")
    }
    # Time Variable checks
    if (length(time_variable) > 1) {
      stop("You can only specify one Time variable")
    }
    if (!(time_variable %in% variables)) {
      stop("Time variable is not part of the provided variables")
    }
    
    # generalized_mode checks
    if(generalized_mode != TRUE){
      stop("generalized_mode should be TRUE for simulating survival data")
    }
    # generalized_mode_model_event/generalized_mode_model_no_event checks
    
    
    data <- data[, variables, drop = FALSE]
    
    bin_variables_surv <- bin_variables[-which(bin_variables == event_variable)]
    
    
    
    # Separating data set to event and no even data set
    dataset_no_event <- data[which(data[,event_variable] == 0),
                             which(names(data) != event_variable)]
    
    dataset_event <- data[which(data[,event_variable] == 1),
                      which(names(data) != event_variable)]
    
    n_samples_no_event <- round(nrow(dataset_no_event)*(n_samples/nrow(data)))
    n_samples_event <- n_samples - n_samples_no_event
    
    # Setting up the arguments to be used in the 2 modgo runs
    .args[["event_variable"]] <- NULL
    .args[["time_variable"]] <- NULL
    .args[["surv_method"]] <- NULL
    
    .args[["generalized_mode_model_no_event"]] <- NULL
    .args[["generalized_mode_model_event"]] <- NULL
    
    .args_corr <- .args
    .args$bin_variables <- bin_variables_surv
    
    
    .args_corr$stop_sim <- TRUE
    .args_corr$generalized_mode_model <- NULL

    correlation_matrix <- do.call(modgo, .args_corr)$covariance_matrix
    # Survival Method number 2
    
    .args$generalized_mode <- TRUE
    # Setting up the arguments to be used for each modgo run
    .args_no_event <- .args_event <- .args
    if(length(generalized_mode_model[,"Variables"])!=0) {
      if(time_variable %in% generalized_mode_model[,"Variables"]){
        stop(paste0("Time variable generalized model should be specified only through generalized_mode_model_no_event and generalized_mode_model_event argument"))
      }
    }
    .args_no_event$data <- dataset_no_event
    .args_no_event$n_samples <- n_samples_no_event
    
    
    if (length(.args_no_event$generalized_mode_model) == 0) {
      
      .args_no_event$generalized_mode_model <- cbind(Variables = time_variable,
                                                     Model = generalized_mode_model_no_event)
      
    } else if (length(.args_no_event$generalized_mode_model) > 0){
      .args_no_event$generalized_mode_model <- rbind(generalized_mode_model,
                                                     c(time_variable, generalized_mode_model_no_event))
    }
    .args_event$data <- dataset_event
    .args_event$n_samples <- n_samples_event
    
    
    if (length(.args_event$generalized_mode_model) == 0) {
      .args_event$generalized_mode_model <- cbind(Variables = time_variable,
                                                  Model = generalized_mode_model_event)
    } else if(length(.args_event$generalized_mode_model) > 0){
      .args_event$generalized_mode_model <- rbind(generalized_mode_model,
                                                     c(time_variable, generalized_mode_model_event))
    }
    # Modgo runs and calculate correlation for all variables
    
    add_status <- function(x, status) {
      
      res <- cbind(x, status)
      colnames(res)[which(colnames(res) == "status")] <- event_variable
      return(res)
    }
    
    if(surv_method == 2){
      
      .args_no_event_cov_matrix <- .args_no_event
      .args_event_cov_matrix <- .args_event
      
      .args_no_event_cov_matrix$stop_sim <- 
        .args_event_cov_matrix$stop_sim <- 
        TRUE
      
      .args_no_event$sigma <- .args_event$sigma <- correlation_matrix[rownames(correlation_matrix)!= event_variable,
                                                                      colnames(correlation_matrix)!= event_variable]
      
      .args_no_event$sigma[,time_variable] <-
        .args_no_event$sigma[time_variable,] <- 
        do.call(modgo,
                .args_no_event_cov_matrix)$covariance_matrix[time_variable,]
      
      .args_event$sigma[,time_variable] <-
        .args_event$sigma[time_variable,] <- 
        do.call(modgo,
                .args_event_cov_matrix)$covariance_matrix[time_variable,]
      
    }
    modgo_no_event <- do.call(modgo, .args_no_event)
    status <- rep(0,dim(modgo_no_event$simulated_data[[1]])[1])
    modgo_no_event$simulated_data <- lapply(modgo_no_event$simulated_data,
                                            FUN = function(x) add_status(x, status))
    

    modgo_event <- do.call(modgo, .args_event)
    status <- rep(1,dim(modgo_event$simulated_data[[1]])[1])
    modgo_event$simulated_data <- lapply(modgo_event$simulated_data,
                                         FUN = function(x) add_status(x, status))
    
    Simulated_Datasets <- lapply(1:nrep,
                                 FUN = function(x) x = rbind(modgo_no_event$simulated_data[[x]],
                                                             modgo_event$simulated_data[[x]]))
    Correlations_matrices <- lapply(Simulated_Datasets,
                                 FUN = function(x) x = cor(x))
    Correlations_matrices[["Mean"]] <- Reduce('+', Correlations_matrices)/nrep
    
    results <- modgo_event
    
    results[["original_data"]] <- data
    results[["simulated_data"]] <- Simulated_Datasets
    results[["correlations"]] <- Correlations_matrices
    results[["correlations_no_event"]] <- modgo_no_event$correlations
    results[["correlations_event"]] <- modgo_event$correlations
    
    results[["covariance_matrix"]] <- c()
    results[["covariance_matrix"]]$no_event <- modgo_no_event$covariance_matrix
    results[["covariance_matrix"]]$event <- modgo_event$covariance_matrix
    
    results[["presim_sigma"]] <- correlation_matrix
    results[["presim_sigma_no_event"]] <- modgo_no_event$presim_sigma
    results[["presim_sigma_event"]] <- modgo_event$presim_sigma
    
    results[["generalized_mode_lmbds"]] <- c()
    results[["generalized_mode_lmbds"]]$no_event <- modgo_no_event$generalized_mode_lmbds
    results[["generalized_mode_lmbds"]]$event <- modgo_event$generalized_mode_lmbds
    
    results[["bin_variables"]] <- c(.args$bin_variables, event_variable)
    results[["event_variable"]] <- event_variable
    results[["time_variable"]] <- time_variable
    
    return(results)
  }
    