#' Calculate Sigma with the help of polychoric and polyserial functions
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
#' @return A correlation matrix.
#' @author Francisco M. Ojeda, George Koliopanos
#' @keywords Normal rank transformation



Sigma_calculation <- function (data,
                               variables,
                               bin_variables,
                               categ_variables,
                               ties_method) {
  
  df_rbi <-
    as.data.frame(do.call(cbind, lapply(data, function (x)
      rbi_normal_transform(x, ties_method = ties_method))))
  
  Sigma <- cov(df_rbi)
  diag(Sigma) <- 1
  
  # Check if categorical variable has
  # more than 7 values. If so it will be treated as continuous
  categ_var <- NULL
  for (ct_var in categ_variables) {
    if (length(unique(data[, ct_var])) > 7) {
      message(
        paste0(
          ct_var,
          " will be treated as continuous due to having more than 7 unique values."
        )
      )
    } else {
      categ_var <- c(categ_var, ct_var)
    }
  }
  
  subst_list <-
    vector(mode = "list", length = length(categ_var))
  names(subst_list) <- categ_var
  # Change values of categorical variables to be consecutive
  # This step is needed for the calculation of polychoric and polyserial
  for (cate in categ_var) {
    unq_values <- sort(unique(data[, cate]))
    
    values = unq_values
    subst_list[[cate]] <- values
    fake_data <- vector(length = length(data[, cate]))
    
    for (i in c(1:length(unq_values))) {
      fake_data[which(data[, cate] == values[i])] <- i
      
    }
    
    data[, cate] <-  fake_data
  }
  oldw <- getOption("warn")
  
  options(warn = -1)
  
  
  Sigma <- suppressWarnings(suppressMessages(Sigma_transformation(
    data = data,
    data_z = df_rbi,
    Sigma,
    variables = variables,
    bin_variables = bin_variables,
    categ_variables = categ_var
  )))
  
  options(warn = oldw)
  
  
  # Transform categorical variables to their original values
  for (cate in categ_var) {
    values <- subst_list[[cate]]
    fake_data <- vector(length = length(data[, cate]))
    for (i in c(1:length(values))) {
      fake_data[which(data[, cate] == i)] <- values[i]
      
    }
    data[, cate] <-  fake_data
  }
  return(Sigma)
}