#' Check Arguments
#'
#' Check that the arguments are following
#' the corresponding conditions
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
#' @param generalized_mode A logical value indicating if you want to use generalized 
#' distribution to simulate your data
#' @param generalized_mode_model A matrix that contains two columns named "Variable" and
#' "Model". This matrix can be used only if a generalized_mode_model argument is
#' provided. It specifies what model should be used for each Variable.
#' Model values should be "RMFMKL", "RPRS", "STAR" or a combination of them,
#' e.g. "RMFMKL-RPRS" or "STAR-STAR", in case the use wants a bimodal simulation.
#' The user can select Generalised Poisson model for poisson variables,
#' but this model cannot be included in bimodal simulation.
#' @param generalized_mode_lmbds A matrix that contains lmbds values for each of the
#' variables of the data set to be used for either Generalized Lambda Distribution
#' Generalized Poisson Distribution or setting up thresholds
#' @param new_mean_sd A matrix that contains two columns named
#' "Mean" and "SD" that the user specifies desired Means and Standard Deviations
#' in the simulated data sets for specific continues variables. The variables
#' must be declared as ROWNAMES in the matrix
#' @author Francisco M. Ojeda, George Koliopanos
#' @export

checkArguments <-
  function(data=NULL,
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
    # Include only selected variables in the original dataset
    data <- data[, variables, drop = FALSE]
    # Find the continuous variables
    continuous_var <-
      setdiff(variables, c(bin_variables, categ_variables))
    
    # Set GLD/GPD accepted methods
    methods <- c("rmfmkl", "rprs", "star")
    methods <- c(methods,
                 apply(expand.grid(methods, methods), 1, paste, collapse =
                         "-"))
    methods <- c(methods, "gp")
    
    #Check input
    if (is.null(data)){
      message("Data set is not provided")
      if(is.null(sigma)){
        stop("Since Data set is not provided, you need to provide the Correlation Matrix(sigma)")
      }
      if(generalized_mode != TRUE){
        stop("Since Data set is not provided, you need to provide generalized_mode should be TRUE")
      }
      if(is.null(generalized_mode_lmbds)){
        stop("Since Data set is not provided, you need to provide generalized_mode_lmbds")
      }
      if(is.null(n_samples)){
        stop("Since Data set is not provided, you need to provide samples size")
      }
    }
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
    # Check that count variables are part of the catgorical variables
    if (!all(count_variables %in% categ_variables)) {
      stop("Count variables are not part of the categorical variables")
    }
    # Check that count variables are used only when generlized mode is true
    if (!is.null(count_variables) & generalized_mode != TRUE) {
      stop("Count variables can only be simulated with as Poisson variables when generalized_mode = TRUE")
    }
    # Check that data is a data frame if dataset is offereed
    if (!is.null(data)){
      if (!is.data.frame(data)) {
        stop("Data is not a data frame")
      }
      # Check that data does not contain NAs
      if (any(is.na(data))) {
        stop("Data set contains NA's")
      }
      # Check that data has at least 2 columns
      if (dim(data)[2] < 2) {
        stop("Data set must contain at least 2 columns")
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
    
    # Check that generalized_mode is either TRUE or FALSE
    if (!(is.logical(generalized_mode))){
      stop("generalized_mode should be a TRUE or FALSE")
    }
    # Check that generalized_mode_model is a matrix
    if ((!is.matrix(generalized_mode_model) &
         length(generalized_mode_model) != 0)) {
      stop("generalized_mode_model should be a matrix")
    }
    # Check that the column names of generalized_mode_model are Variables and Model
    if ((!all(colnames(generalized_mode_model) %in% c("Variables", "Model"))
         & length(generalized_mode_model) != 0)) {
      stop("generalized_mode_model colnames should be Variables and Model")
    }
    # Check that Variables from generalized_mode_model are part of the 
    # declared variables
    if (!all(generalized_mode_model[, "Variables"] %in% variables)) {
      stop("generalized_mode_model Variables should be a part of provided variables")
    }
    # Check that Model from generalized_mode_model are part of the available methods
    if (!all(generalized_mode_model[, "Model"] %in% methods)) {
      stop(
        paste0("All generalized_mode_model Model should be part of "),
        paste(methods, collapse = ", ")
      )
    }
    # Check that generalized_mode_model does not contain duplicate variables or 
    # blank values
    if(length(generalized_mode_model) != 0){
      if (length(unique(generalized_mode_model[, "Variables"])) != dim(generalized_mode_model)[1]) {
        stop("generalized_mode_model Variables contains either duplicate variable or blank value")
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
      stop("Pertrubation variables are not part of the provided
         continuous variables")
    }
    # Check that var_infl includes only continuous variables
    if (length(var_infl) > 0 &
        !all(names(var_infl) %in% continuous_var)) {
      stop("Variance inflation variables are not part of the provided
         continuous variables")
    }
    # Check that infl_cov_stable is TRUE or FALSE
    if (!is.logical(infl_cov_stable)) {
      stop("Inflation variance stable should be TRUE/FALSE")
    }
    # Check that the user select either one of prtr_vec or var_infl and not both
    if (any(names(pertr_vec) %in% names(var_infl))) {
      stop("Perturb vector cannot have common variables with variance inflation vector")
    }
    # Ceck that the change amount and the pair of change covariance has been
    # both set
    if (change_amount != 0 & length(change_cov) == 0) {
      stop(
        "You need to provide a pair of variables(change_cov)
          to change their covariance values"
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

