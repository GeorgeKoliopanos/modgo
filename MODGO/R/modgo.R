#' MOck Data GeneratiOn
#' 
#' \code{modgo} Create mock dataset from a real one by using 
#' ranked inverse normal transformation.
#' 
#' [We can add more details here, once the paper is almost finished]
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
#' @param bin_variables  a character vector listing the binary variables.
#' @param categ_variables a character vector listing the ordinal categorical 
#' variables.
#' @param nrep number of repetitions.
#' @param noise_mu Logical value if you want to apply noise to  
#' multivariate mean. Default: FALSE
#' @param pertr_vec A vector of variables that the user wants to perturb
#' @param nprod Number of rows of each simulated dataset. Default is
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
#' @param thresh_right A numeric vector indicating the lower limits of the 
#' the specified variables
#' @param thresh_left A numeric vector indicating the upper limits of the 
#' the specified variables
#' @param thresh_force A logical value indicating if you want to force threshold
#' in case the proportion of samples that can surpass the threshold are less 
#' than 10\%
#' @param var_prop A named vector that provides a  proportion of 
#'  value=1 for a specific binary variable(=name of the vector) that will be
#'  the proportion of this value in the simulated data sets.[this may increase
#'  execution time drastically]
#' @return A list with the following components:
#' \item{SimulatedData}{A list of data frames containing the simulated data.}
#' \item{OriginalData}{A data frame with the input data.}
#' \item{Correlations}{a list of correlation matrices. The ith element is the
#' correlation matrix for the ith simulated dataset. The \code{(repn + 1)}th
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
#' @examples 
#' data("Cleveland",package="modgo")
#' test_modgo <- modgo(data = Cleveland,
#'      bin_variables = c("CAD","HighFastBloodSugar","Sex","ExInducedAngina"),
#'      categ_variables =c("Chestpaintype"))
#' @export
#' @importFrom Matrix nearPD

modgo <- function(data,ties_method=  "max", variables= colnames(data),
                  bin_variables= c(),categ_variables= c(),
                  nprod= nrow(data),sigma= c(),nrep= 100,
                  noise_mu= FALSE, pertr_vec= c(),
                  change_cov= c(),change_amount= 0,seed= 1,
                  thresh_var= c(), thresh_force = FALSE, 
                  var_prop= c()) {
  
  if (!is.na(seed)){
  # Setting Seed
  set.seed(seed)
  
  }

  #Check input  
  
  if (!is.data.frame(data)){
    stop("Data is not a data frame")
  }
  
  
  if (!all(variables %in% colnames(data))){
  
  stop("Variables are not column names of data ")
  
  }
  
  if (length(sigma)>0 & (colnames(sigma)!= variables || 
                         rownames(sigma)!= variables)){
    
    stop("Sigma column and row names and not exactly equal to variable names")
    
  }
  # Include only selected variables in the original dataset
  data <- data[,variables]
  
  if(!all(apply(data, 2, function(x) is.numeric(x)))){
    
    stop("Data should only contain numerical values")
    
  }
  if (!all(bin_variables %in% variables)){
  
  stop("Binary variables are not part of the provided variables")
  
  }
  if (!all(thresh_var[,1] %in% variables)){
    
    stop("Threshold variables are not part of the provided variables")
    
  }
  
  # Find the continuous variables
  continuous_var <- setdiff(variables,c(bin_variables,categ_variables))
  
  if (length(pertr_vec) > 0 && !all(names(pertr_vec) %in% continuous_var) ){
    
    stop("Pertrubation variables are not part of the provided 
         continuous variables")
    
  }
  
  if (change_amount!=0 && length(change_cov)==0){
    
  stop("You need to provide a pair of variables(change_cov) 
       to change their covariance values")
  }
  if (change_amount==0 && length(change_cov)==2){
    
    stop("You need to provide an amount of change (change_amount)")
  }
  
  if (!(length(change_cov) %in% c(0,2))) {
    
    stop("You need to provide two variables")
  }
  # Noise in multivariate distributions centers
  if (noise_mu == TRUE){
    
    ns <- rnorm(ncol(data),mean = 0,sd=1)
    
    
  } else {
    ns <- matrix(0,nrow = 1,ncol = ncol(data))
    
  }
  
  if (length(var_prop)>0 & !all(names(var_prop) %in% bin_variables)){
    stop("Name of var_prop should be a binary variable")
    
  }
  if (length(var_prop)>0 & (var_prop>1 || var_prop <0)){
    stop("Variable proportion should be between 0 and 1")
    
  }
  
  if (length(var_prop)>1 ){
    stop("You cannot set proportion for more than 1 variable")
    
  }
  if (length(var_prop) ==1 & length(thresh_var[,1]) >0 ){
    stop("You cannot set a variable proportion and 
         a variable threshold at the same run ")
    
  }
  OriginalData <- data
  
  # Normal run in case user does not provide a sigma
  if(length(sigma)==0){
    
  # Rank transofrmation of each column of the data set
  df_rbi <- data
  for(j in 1:ncol(df_rbi)) {
    df_rbi[[j]] <- rbi_normal_transform(data[[j]],ties_method = ties_method)
  }
  Sigma <- cov(df_rbi)
  diag(Sigma) <- 1
  
  categ_var <- c()
  for (ct_var in categ_variables){
    
    if(length(unique(data[,ct_var])) > 7){ 
      print(
      paste0(ct_var,
     " will be treated as continuous due to having more than 7 unique values."))
      }else {
                     categ_var <- c(categ_var,ct_var)
  }
  }
    
  subst_list <- vector(mode="list",length = length(categ_var))
  names(subst_list) <- categ_var
  for (cate in categ_var){
    
    unq_values <- sort(unique(data[,cate]))

    values = unq_values
    subst_list[[cate]] <- values
    
    
    for (i in c(1:length(unq_values))){
      
      data[which(data[,cate]==values[i]),cate] <- i
      
    } 
    
  }
  
  oldw <- getOption("warn")
  
  options(warn = -1)
  
  
  Sigma <- suppressMessages(Sigma_transformation(
    data=data,data_z = df_rbi,Sigma,
    variables = variables,
    bin_variables = bin_variables,
    categ_variables=categ_var))
 
  options(warn = oldw)
 
  
  
  for (cate in categ_var){
    
    values <- subst_list[[cate]]
    for (i in c(length(values):1)) {
      
      data[which(data[,cate]==i),cate] <- values[i]
      
      
      
    }
    
  }
  }else {Sigma <- sigma}
  
# Change specified pairs of covariance matrix by specified amount
 if (change_amount!=0 && length(change_cov)==2) {
   Sigma[change_cov[1],change_cov[2]] <-  
     Sigma[change_cov[2],change_cov[1]] <- 
     Sigma[change_cov[2],change_cov[1]] + change_amount 
   
 }
  
  
  ## Thresholds
  if (length(thresh_var[,1]) >0){
    mt_sim <- MASS::mvrnorm(n = nprod, 
                            mu = rep(0, ncol(data)) +ns,
                            Sigma = Sigma)
    df_sim <- data.frame(mt_sim)
    names(df_sim) <- names(data)
    #Inverse transformation of each variable
    for(j in 1:ncol(df_sim)) {
      df_sim[[j]] <- rbi_normal_transform_inv(df_sim[[j]], data[[j]])
      
    }
    
    for (i in c(1:length(thresh_var[,1]))){
      if(is.na(thresh_var[i,2]))
      {low_thresh <-
          min(df_sim[thresh_var[i,1]])-1}else{low_thresh <-thresh_var[i,2]}
      if(is.na(thresh_var[i,3]))
      {up_thresh <-
          max(df_sim[thresh_var[i,1]])+1}else{up_thresh <-thresh_var[i,3]}
      
      
      df_sim <- df_sim[which(df_sim[thresh_var[i,1]]< up_thresh &
                               df_sim[thresh_var[i,1]]> low_thresh),]
    }
    
      if(length(df_sim[,1])/nprod <0.1 & thresh_force == FALSE ){
      stop("The propotion of simulated samples passing
           the threshold is less than 10% ")
      }else if(length(df_sim[,1])==0){
        stop("The propotion of simulated samples passing
           the threshold is 0% ")
      }else{thresh_multi <- 1.1*nprod/length(df_sim[,1])}
    
    
   }else{thresh_multi=1}
  
  # Binary variables proportions
  if (length(var_prop) ==1 ){
    
    mt_sim <- MASS::mvrnorm(n = nprod, 
                            mu = rep(0, ncol(data)) +ns,
                            Sigma = Sigma)
    df_sim <- data.frame(mt_sim)
    names(df_sim) <- names(data)
    #Inverse transformation of each variable
    for(j in 1:ncol(df_sim)) {
      df_sim[[j]] <- rbi_normal_transform_inv(df_sim[[j]], data[[j]])
      
    }
    
    counts_1 <- length(which(df_sim[,names(var_prop)]==1))
    counts_0 <- length(which(df_sim[,names(var_prop)]==0))
    req_1 <- var_prop * nprod
    req_0 <- (1-var_prop) * nprod
    
    if (req_1 >= counts_1){
    thresh_multi <-  (1.1*req_1)/counts_1
      
    }else if (req_1 < counts_1){
      thresh_multi <-  (1.1*req_0)/counts_0
      
    }
    
    
  }
  
  # Check Sigma if it is positive definite
  if(!all(eigen(Sigma)$value > 0)){
    
    print("Covariance matrix is not positive definite.") 
    print("It will be replaced with nearest positive definite matrix")
    Sigma <- Matrix::nearPD(Sigma,corr = TRUE)$mat
  }
  
  # Starting loop for many repetions
  
  Correlations <- vector(mode = "list", length = nrep+1)
  mean_corr <- matrix(0,nrow =ncol(data),ncol = ncol(data))
  
  SimulatedData <- vector(mode = "list", length = nrep)
  i <- 1
  counter <- 0
  # Loop for creating new datasets and obtaining their mean correlations
  while (i < nrep+1){
    
  mt_sim <- MASS::mvrnorm(n = ceiling(nprod*thresh_multi), 
                          mu = rep(0, ncol(data)) + ns,
                          Sigma = Sigma)
  
  
  df_sim <- data.frame(mt_sim)
  names(df_sim) <- names(data)
  #Inverse transformation of each variable
  for(j in 1:ncol(df_sim)) {
    df_sim[[j]] <- rbi_normal_transform_inv(df_sim[[j]], data[[j]])
    if(colnames(data)[[j]] %in% names(pertr_vec)){
      
      p <- pertr_vec[which(
        names(pertr_vec)==colnames(data)[[j]])]
      
      df_sim[[j]] <- (df_sim[[j]]*sqrt(1-p)) + 
        rnorm(length(df_sim[[j]]),mean = 0,sd=sd(df_sim[[j]])*sqrt(p))
      
    }
  }
  #Threshold process
  if(length(thresh_var[,1])>0){
  for (j in c(1:length(thresh_var[,1]))){
    if(is.na(thresh_var[j,2]))
    {low_thresh <-
      min(df_sim[thresh_var[j,1]])-1}else{low_thresh <-thresh_var[j,2]}
    if(is.na(thresh_var[j,3]))
    {up_thresh <-
      max(df_sim[thresh_var[j,1]])+1}else{up_thresh <-thresh_var[j,3]}
    
    
    df_sim <- df_sim[which(df_sim[thresh_var[j,1]]< up_thresh &
                             df_sim[thresh_var[j,1]]> low_thresh),]
  }
  
  
  }
  #Proportion process
  if(length(var_prop)==1){

      rounded_length <- round(nprod*var_prop)
      df_sim_1 <- df_sim[which(df_sim[,names(var_prop)] == 1),]
      df_sim_0 <- df_sim[which(df_sim[,names(var_prop)] == 0),]
      df_sim <- rbind(df_sim_1[c(1:rounded_length),],
                      df_sim_0[c(1:(nprod-rounded_length)),])
      
      
      
      
    
  }
  
  if(length(df_sim[,1]) < nprod || is.null(df_sim)){
    i <- i-1
    counter <-counter +1
  } else {
    df_sim <- df_sim[c(1:nprod),]
    
    #Correlation calculation
    Correlations[[i]] <- cor(df_sim)
    SimulatedData[[i]] <- df_sim
    
    #Mean correlation calculation
    mean_corr <- mean_corr + (cor(df_sim)/nrep)
  }
  
  if(counter > 10*nrep){
    stop("Could not create simulated datasets with this specific thresholds")
  } 
  
  i <- i + 1
  }
  
  
  Correlations[[nrep+1]] <- mean_corr
  names(Correlations) <- c(paste0("rep",seq(1:nrep)),"Mean")
  
  results <-list(SimulatedData,OriginalData,Correlations,
                 bin_variables,categ_variables,Sigma,seed,nprod,nrep)
  names(results) <-c("SimulatedData","OriginalData","Correlations",
                     "Binary_variables","Categorical_variables",
                     "Covariance_Matrix","Seed","Samples_Produced",
                     "Sim_Dataset_Number")
  return(results)
}

