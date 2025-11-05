#' Check Arguments
#'
#' This function is used internally by \code{modgo_svy} to check the correctness 
#' of the arguments passed to it.  
#' 
#' All variables passed to \code{modgo_svy} should be of class 
#' double or integer. This includes the variables passed to the parameter
#' \code{categ_variables}. The character vector \code{variables}, indicating
#' the variables in \code{data} to be used in the simulation, should 
#' contain at least two variables. The variables in \code{variables} not present
#' in \code{bin_variables} nor \code{categ_variables} will be treated as 
#' continuous variables.
#' 
#' @inheritParams modgo_svy
#' @author Francisco M. Ojeda, George Koliopanos
#' @export

checkArguments_svy <-
  function(design = NULL,
           variables = NULL,
           bin_variables = NULL,
           categ_variables = NULL,
           n_samples = nrow(design[["variables"]]),
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
           multi_sugg_prop = NULL) {
    
    # Find the continuous variables
    continuous_var <-
      setdiff(variables, c(bin_variables, categ_variables))
    
    #Check input
    if(is.null(design)) stop("design needs to be provided")
    
    # Check that number of repetitions is an Integer number above 0
    if (nrep < 1) {
      stop("nrep should be larger than 0")
    }
    if (!nrep == as.integer(nrep)) {
      stop("nrep should be an integer Number")
    }
    
    # Check that number of samples is an Integer number above 0
    if (n_samples < 1) {
      stop("n_samples should be larger than 0")
    }
    
    if (!n_samples == as.integer(n_samples)) {
      stop("n_samples should be an integer Number")
    }
    
    # Check that binary variables are part of the declared variables
    if (!all(bin_variables %in% variables)) {
      stop("Binary variables are not part of the provided variables")
    }
    
    # Check that categorical variables are part of the declared variables
    if (!all(categ_variables %in% variables)) {
      stop("Categorical variables are not part of the provided variables")
    }
    
    # Check that data is a data frame if dataset is offereed
    if (!is.null(design)){
      if (!"survey.design2" %in% class(design)) {
        stop("design is not of class 'survey.design2'")
      }
      
      # Check that data does not contain NAs
      data <- design[["variables"]][, variables]
      if (any(is.na(data))) {
        stop("Dataset contains NA's")
      }
      # Check that data has at least 2 columns
      if (dim(data)[2] < 2) {
        stop("Dataset must contain at least 2 columns")
      }
      # Check that declared variables are part of the column names of data
      if (!all(variables %in% colnames(data))) {
        stop("Variables are not column names of data ")
      }
      # Check that all binary variables have only 0 and 1 as values & that they 
      # have both values present
      for (bin in bin_variables) {
        if (length(which(data[, bin] %in% c(0, 1))) < length(data[, bin])) {
          stop(paste0("Binary variable ", bin, " is neither 0 nor 1"))
        }
        if (length(unique(data[, bin])) == 1) {
          stop(paste0("Binary variable ", bin, " has only value"))
        }
      }
      # Check that data contains only numerical values
      if (!all(apply(data, 2, function(x)
        is.numeric(x)))) {
        stop("Data should only contain numerical values")
      }
    }
    
    # Check that Sigma is not a 0 length matrix
    if (length(sigma) > 0){
      if(all(colnames(sigma) != variables) | all(rownames(sigma) != variables)) {
        stop("Sigma column and row names and not exactly equal to variable names")
      }
    }
    
    
    # Check that the row names of new_mean_sd are part of the declared variables
    if (!all(rownames(new_mean_sd) %in% variables)) {
      stop("Rownames of new_mean_sd are not part of the provided variable")
    }
    # Check that new_mean_sd column names are Mean and SD
    if (length(new_mean_sd) != 0) {
      if (!all(colnames(new_mean_sd) %in% c("Mean", "SD")) |
          dim(new_mean_sd)[2] != 2) {
        stop("Colnames of new_mean_sd are not Mean and SD")
      }
    }
    # Check that multi_sugg_prop are part of the declared binary variables
    if (!all(names(multi_sugg_prop) %in% bin_variables)) {
      stop("Names of multi_sugg_prop are not part of the provided binary variables")
    }
    # Check that thresh_var are part of the declared variables
    if (!all(thresh_var[, 1] %in% variables)) {
      stop("Threshold variables are not part of the provided variables")
    }
    # Check pertr_vec includes only continuous variables
    if (length(pertr_vec) > 0 &
        !all(names(pertr_vec) %in% continuous_var)) {
      stop(paste0("Perturbation variables are not part of the provided ",
                  "continuous variables"))
    }
    # Check that var_infl includes only continuous variables
    if (length(var_infl) > 0 &
        !all(names(var_infl) %in% continuous_var)) {
      stop(paste0("Variance inflation variables are not part of the provided ",
                  "continuous variables"))
    }
    # Check that infl_cov_stable is TRUE or FALSE
    if (!is.logical(infl_cov_stable)) {
      stop("Inflation variance stable should be TRUE/FALSE")
    }
    # Check that the user select either one of prtr_vec or var_infl and not both
    if (any(names(pertr_vec) %in% names(var_infl))) {
      stop(paste0("Perturbation vector cannot have common variables with ",
                  "variance inflation vector"))
    }
    # Check that the change amount and the pair of change covariance has been
    # both set
    if (change_amount != 0 & length(change_cov) == 0) {
      stop(
        paste0("You need to provide a pair of variables(change_cov) to ",
               "change their covariance values")
      )
    }
    # Check if the user has specify the amount for a specific variable he wants
    #  to change
    if (change_amount == 0 & length(change_cov) == 2) {
      stop("You need to provide an amount of change (change_amount)")
    }
    # Check if for change covariance arguments the user has define a pair of 
    # variables
    if (!(length(change_cov) %in% c(0, 2))) {
      stop("You need to provide two variables")
    }
    # Check if the variable that it is used for a set up proportion is part
    # of the binary variables
    if (length(var_prop) > 0 &
        !all(names(var_prop) %in% bin_variables)) {
      stop("Name of var_prop should be a binary variable")
    }
    # Check if variable proportion is between the allowed limit 
    if (length(var_prop) > 0 ){ 
      if(var_prop > 1 | var_prop < 0){
        stop("Variable proportion should be between 0 and 1")
      }
    }
    # Check if user has set a proportion for more than 1 variables
    if (length(var_prop) > 1) {
      stop("You cannot set proportion for more than 1 variable")
    }
    # Check if user has declare a binary variable proportion and
    # a continuous variable threshold at the same time
    if (length(var_prop) == 1 & length(thresh_var[, 1]) > 0) {
      stop("You cannot set a variable proportion and
         a variable threshold at the same run ")
    }
    
    # Check Sigma if it contains NA's
    if (any(is.na(sigma))) {
      stop("Correlation of data contains NA's")
    }
}

