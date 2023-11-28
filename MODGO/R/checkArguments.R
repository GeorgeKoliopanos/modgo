#' Check Arguments
#'
#' Check that the arguments are following
#' the corresponding conditions
#' #' @param data a data frame containing the data whose characteristics are to be
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
#' @param gener_var A vector indicating the which variables should
#' be simulated using generalized lambda distribution
#' @param gener_var_model A matrix that contains two columns named "Variable" and
#' "Model". This matrix can be used only if a gener_var_model argument is
#' provided. It specifies what model should be used for each Variable.
#' Model values should be "RMFMKL", "RPRS", "STAR" or a combination of them,
#' e.g. "RMFMKL-RPRS" or "STAR-STAR", in case the use wants a bimodal simulation.
#' The user can select Generalised Poisson model for poisson variabes,
#' but this model cannot be included in bimodal simulation.
#' @param new_mean_sd A matrix that contains two columns named
#' "Mean" and "SD" that the user specifies desired Means and Standard Deviations
#' in the simulated data sets for specific continues variables. The variables
#' must be declared as ROWNAMES in the matrix
#' @export

checkArguments <-
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
           gener_var = NULL,
           gener_var_model = NULL,
           gener_var_lmbds = NULL) {
    #.args <- as.list(match.call()[-1])
    
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
    if (!is.data.frame(data)) {
      stop("Data is not a data frame")
    }
    if (any(is.na(data))) {
      stop("Data set contains NA's")
    }
    if (dim(data)[2] < 2) {
      stop("Data set must contain at least 2 columns")
    }
    
    if (!all(variables %in% colnames(data))) {
      stop("Variables are not column names of data ")
      
    }
    if (!all(bin_variables %in% variables)) {
      stop("Binary variables are not part of the provided variables")
      
    }
    
    for (bin in bin_variables) {
      if (length(which(data[, bin] %in% c(0, 1))) < length(data[, bin])) {
        stop(paste0("Binary variable ", bin, " is neither 0 nor 1"))
        
      }
      if (length(unique(data[, bin])) == 1) {
        stop(paste0("Binary variable ", bin, " has only value"))
        
      }
    }
    if (!all(categ_variables %in% variables)) {
      stop("Categorical variables are not part of the provided variables")
      
    }
    if (!all(count_variables %in% categ_variables)) {
      stop("Count variables are not part of the categorical variables")
      
    }
    
    
    if (!(is.vector(gener_var) | length(gener_var) == 0)) {
      stop("gener_var should be a vector")
    }
    
    if (!all(gener_var %in% c(continuous_var, categ_variables))) {
      stop("gener_var should be either in categorical or continuous variables")
    }
    if ((!is.matrix(gener_var_model) &
         length(gener_var_model) != 0)) {
      stop("gener_var_model should be a matrix")
    }
    if ((!all(colnames(gener_var_model) %in% c("Variables", "Model"))
         & length(gener_var_model) != 0)) {
      stop("gener_var_model colnames should be Variables and Model")
    }
    if (!all(gener_var_model[, "Variables"] %in% gener_var)) {
      stop("gener_var_model Variables should be a part of gener_var")
    }
    if (!all(gener_var_model[, "Model"] %in% methods)) {
      stop(
        paste0("All gener_var_model Model should be part of "),
        paste(methods, collapse = ", ")
      )
    }
    if(length(gener_var_model) != 0){
      if (length(unique(gener_var_model[, "Variables"])) != dim(gener_var_model)[1]) {
        stop("gener_var_model Variables contains either duplicate variable or blank value")
      }
    }
    
    if (length(sigma) > 0){
      if(all(colnames(sigma) != variables) | all(rownames(sigma) != variables)) {
        stop("Sigma column and row names and not exactly equal to variable names")
      
      }
    }
    # Include only selected variables in the original dataset
    data <- data[, variables, drop = FALSE]
    
    if (!all(apply(data, 2, function(x)
      is.numeric(x)))) {
      stop("Data should only contain numerical values")
      
    }
    
    if (!all(rownames(new_mean_sd) %in% variables)) {
      stop("Rownames of new_mean_sd are not part of the provided variable")
      
    }
    if (length(new_mean_sd) != 0) {
      if (!all(colnames(new_mean_sd) %in% c("Mean", "SD")) |
          dim(new_mean_sd)[2] != 2) {
        stop("Colnames of new_mean_sd are not Mean and SD")
        
      }
    }
    if (!all(names(multi_sugg_prop) %in% bin_variables)) {
      stop("Names of multi_sugg_prop are not part of the provided binary variables")
      
    }
    if (!all(thresh_var[, 1] %in% variables)) {
      stop("Threshold variables are not part of the provided variables")
      
    }
    
    
    
    if (length(pertr_vec) > 0 &
        !all(names(pertr_vec) %in% continuous_var)) {
      stop("Pertrubation variables are not part of the provided
         continuous variables")
      
    }
    
    if (length(var_infl) > 0 &
        !all(names(var_infl) %in% continuous_var)) {
      stop("Variance inflation variables are not part of the provided
         continuous variables")
      
    }
    if (!is.logical(infl_cov_stable)) {
      stop("Inflation variance stable should be TRUE/FALSE")
      
    }
    
    if (any(names(pertr_vec) %in% names(var_infl))) {
      stop("Perturb vector cannot have common variables with variance inflation vector")
      
    }
    
    if (change_amount != 0 & length(change_cov) == 0) {
      stop(
        "You need to provide a pair of variables(change_cov)
          to change their covariance values"
      )
    }
    
    if (change_amount == 0 & length(change_cov) == 2) {
      stop("You need to provide an amount of change (change_amount)")
      
    }
    
    if (!(length(change_cov) %in% c(0, 2))) {
      stop("You need to provide two variables")
      
    }
    
    if (length(var_prop) > 0 &
        !all(names(var_prop) %in% bin_variables)) {
      stop("Name of var_prop should be a binary variable")
      
    }
    
    if (length(var_prop) > 0 ){ 
      if(var_prop > 1 | var_prop < 0){
        stop("Variable proportion should be between 0 and 1")
      }
    }
    
    if (length(var_prop) > 1) {
      stop("You cannot set proportion for more than 1 variable")
      
    }
    
    if (length(var_prop) == 1 & length(thresh_var[, 1]) > 0) {
      stop("You cannot set a variable proportion and
         a variable threshold at the same run ")
      
    }
    
    # Check Sigma if it contains NA's
    if (any(is.na(sigma))) {
      stop("Correlation of data contains NA's")
    }
  }

