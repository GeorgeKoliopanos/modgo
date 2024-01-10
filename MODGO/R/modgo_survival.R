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
#' @param gener_var A logical value indicating if generalized lambda/poisson
#'  distributions or set up thresholds will be used to generate the simulated values
#' @param gener_var_model_no_event A matrix that contains two columns named "Variable" and
#' "Model" and it is to be used for the non-event data set(event = 0). This matrix can be used only if a gener_var_model argument is
#' provided. It specifies what model should be used for each Variable.
#' Model values should be "rmfmkl", "rprs", "star" or a combination of them,
#' e.g. "rmfmkl-rprs" or "star-star", in case the use wants a bimodal simulation.
#' The user can select Generalised Poisson model for poisson variables,
#' but this model cannot be included in bimodal simulation
#' @param gener_var_model_event A matrix that contains two columns named "Variable" and
#' "Model" and it is to be used for the event data set(event = 1). This matrix can be used only if a gener_var_model argument is
#' provided. It specifies what model should be used for each Variable.
#' Model values should be "rmfmkl", "rprs", "star" or a combination of them,
#' e.g. "rmfmkl-rprs" or "star-star", in case the use wants a bimodal simulation.
#' The user can select Generalised Poisson model for poisson variables,
#' but this model cannot be included in bimodal simulation
#' @param gener_var_lmbds A matrix that contains lmbds values for each of the
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
#' \item{SimulatedData}{A list of data frames containing the simulated data.}
#' \item{OriginalData}{A data frame with the input data.}
#' \item{Correlations}{a list of correlation matrices. The ith element is the
#' correlation matrix for the ith simulated dataset. The \code{(repn + 1)}the
#' (last) element of the list is the average of the correlation matrices.}
#' \item{Binary_variables}{character vector listing the binary variables}
#' \item{Categorical_variables}{a character vector listing the ordinal
#' categorical variables}
#' \item{Covariance_Matrix}{Covariance matrix used when generating observations
#' from a multivariate normal distribution.}
#' \item{Seed}{Random seed used.}
#' \item{Samples_Produced}{Number of rows of each simulated dataset.}
#' \item{Sim_Dataset_Number}{Number of simulated datasets produced.}
#' @author Francisco M. Ojeda, George Koliopanos
#' @keywords mock data generation
#' @export
#' @importFrom Matrix nearPD
#' @import dplyr

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
           gener_var = TRUE,
           gener_var_model_event = NULL,
           gener_var_model_no_event = NULL,
           gener_var_lmbds = NULL) {
    
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
    
    # gener_var checks
    if(gener_var != TRUE){
      stop("gener_var should be TRUE for simulating survival data")
    }
    # gener_var_model_event/gener_var_model_no_event checks
    if ((!all(colnames(gener_var_model_no_event) %in% c("Variables", "Model"))
         & length(gener_var_model_no_event) != 0)) {
      stop("gener_var_model_no_event colnames should be Variables and Model")
    }
    if ((!all(colnames(gener_var_model_event) %in% c("Variables", "Model"))
         & length(gener_var_model_event) != 0)) {
      stop("gener_var_model_event colnames should be Variables and Model")
    }
    if(length(gener_var_model_no_event[which(gener_var_model_no_event[,1] == time_variable),2]) == 0){
      gener_var_model_no_event <- rbind(gener_var_model_no_event, c(time_variable, "rprs"))
      colnames(gener_var_model_no_event) <- c("Variables", "Model")
    }
    if(length(gener_var_model_event[which(gener_var_model_event[,1] == time_variable),2]) == 0){
      gener_var_model_event <- rbind(gener_var_model_event, c(time_variable, "rprs"))
      colnames(gener_var_model_event) <- c("Variables", "Model")
    }
    
    
    data <- data[, variables, drop = FALSE]
    
    bin_variables_surv <- bin_variables[-which(bin_variables == event_variable)]
    .args$bin_variables <- bin_variables_surv
    
    
    
    dataset_no_event <- data %>%
      filter(!!as.name(event_variable) == 0) %>%
      select(-!!as.name(event_variable))
    dataset_event <- data %>%
      filter(!!as.name(event_variable) == 1) %>%
      select(-!!as.name(event_variable))
    
    n_samples_no_event <- round(nrow(dataset_no_event)*(n_samples/nrow(data)))
    n_samples_event <- n_samples - n_samples_no_event
    
    # Setting up the arguments to be used in the 2 modgo runs
    .args[["event_variable"]] <- NULL
    .args[["time_variable"]] <- NULL
    .args[["surv_method"]] <- NULL
    
    .args[["gener_var_model_no_event"]] <- NULL
    .args[["gener_var_model_event"]] <- NULL
    
    
    
    .args_corr <- .args
    .args_corr$stop_sim <- TRUE
    correlation_matrix <- do.call(modgo, .args_corr)$Covariance_Matrix
    # Survival Method number 2
    if(surv_method == 2 & length(.args$sigma) == 0){
      .args$sigma <- correlation_matrix[rownames(correlation_matrix)!= event_variable,
                                        colnames(correlation_matrix)!= event_variable]
    }
    
    # Setting up the arguments to be used for each modgo run
    .args_no_event <- .args_event <- .args
    
    .args_no_event$data <- dataset_no_event
    .args_no_event$n_samples <- n_samples_no_event
    .args_no_event$gener_var_model <- gener_var_model_no_event
    
    .args_event$data <- dataset_event
    .args_event$n_samples <- n_samples_event
    .args_event$gener_var_model <- gener_var_model_event
    

    # Calculate correlation for all variables
    modgo_no_event <- do.call(modgo, .args_no_event)
    status <- rep(0,dim(modgo_no_event$SimulatedData[[1]])[1])
    modgo_no_event$SimulatedData <- lapply(modgo_no_event$SimulatedData,
                                           FUN = function(x) x = cbind(x,status))
    
    modgo_event <- do.call(modgo, .args_event)
    status <- rep(1,dim(modgo_event$SimulatedData[[1]])[1])
    modgo_event$SimulatedData <- lapply(modgo_event$SimulatedData,
                                           FUN = function(x) x = cbind(x,status))
    
    Simulated_Datasets <- lapply(1:nrep,
                                 FUN = function(x) x = rbind(modgo_no_event$SimulatedData[[x]],
                                                             modgo_event$SimulatedData[[x]]))
    Correlations_matrices <- lapply(Simulated_Datasets,
                                 FUN = function(x) x = cor(x))
    Correlations_matrices[["Mean"]] <- Reduce('+', Correlations_matrices)/nrep
    
    results <- modgo_event
    
    results[["OriginalData"]] <- data
    results[["SimulatedData"]] <- Simulated_Datasets
    results[["Correlations"]] <- Correlations_matrices
    
    results[["Covariance_Matrix"]] <- c()
    results[["Covariance_Matrix"]]$no_event <- modgo_no_event$Covariance_Matrix
    results[["Covariance_Matrix"]]$event <- modgo_event$Covariance_Matrix
    
    results[["PreSim_Sigma"]] <- correlation_matrix
    results[["PreSim_Sigma_no_event"]] <- modgo_no_event$PreSim_Sigma
    results[["PreSim_Sigma_event"]] <- modgo_event$PreSim_Sigma
    
    results[["Generalized_Lambdas"]] <- c()
    results[["Generalized_Lambdas"]]$no_event <- modgo_no_event$Generalized_Lambdas
    results[["Generalized_Lambdas"]]$event <- modgo_event$Generalized_Lambdas
    
    results[["Binary_variables"]] <- c(.args$bin_variables,event_variable)
    results[["Event_variable"]] <- event_variable
    results[["Time_variable"]] <- time_variable
    
    return(results)
  }
    