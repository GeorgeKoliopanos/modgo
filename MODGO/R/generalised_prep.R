#' Generalized Lambda and Poisson preparation
#' 
#' Prepare the four moments matrix for GLD and GPD
#' @param data a data frame with original variables.
#' @param bin_variables  a character vector listing the binary variables.
#' @param gener_var_model A matrix that contains two columns named "Variables" and
#' "Model". This matrix can be used only if a gener_var_model argument is
#' provided. It specifies what model should be used for each Variable.
#' Model values should be "RMFMKL", "RPRS", "STAR" or a combination of them,
#' e.g. "RMFMKL-RPRS" or "STAR-STAR", in case the use wants a bimodal simulation.
#' The user can select Generalised Poisson model for poisson variabes,
#' but this model cannot be included in bimodal simulation
#' @param multi_sugg_prop A named vector that provides a  proportion of
#'  value=1 for specific binary variables(=name of the vector) that will be
#'  the close to the proportion of this value in the simulated data sets.
#' @return A numeric matrix
#' @author Francisco M. Ojeda, George Koliopanos
#' @export
generalizedMatrix <- function(data,
                              variables = colnames(data),
                              bin_variables = NULL,
                              gener_var_model = NULL,
                              multi_sugg_prop = NULL){
  
  # Prepare gener_var_lmbds matrix
  gener_var_lmbds <- matrix(nrow = 11, ncol = length(variables))
  colnames(gener_var_lmbds) <- variables
  
  for (i in variables){
    if (i %in% gener_var_model[,"Variables"]){
      # If models are provided select the appropriate GLD lambdas calculation
      model <- unlist(strsplit(
        gener_var_model[which(gener_var_model[,"Variables"] == i), 
                        "Model"],
        split = "-"))
      if(length(model) == 2){
        # Bivariate Generalized Lambdas distribution
        biv_lmbds <- tryCatch(GLDEX::fun.auto.bimodal.ml(data[[i]],
                                init1.sel=model[1],
                                init2.sel=model[2],
                                init1=c(-0.25,1.5),
                                init2=c(-0.25,1.5),
                                leap1=3,leap2=3)$par,
                              error = function (e){
                                if(e$message == "non-finite value supplied by optim"){
                                  message(paste0("GLD cannot produce Lambdas with the selection of models: ",paste(model, collapse = "-"), " for variable: ",i))
                                  message(paste0("Default modgo will be used to simulate ", i))
                                  
                                  return(rep(NA, 5))
                                }else{
                                  stop(paste0("Error in the creation of generalised lambdas with bivariate model for variable: ", i))
                                }
                              })
        if (length(biv_lmbds) == 1){
          stop(biv_lmbds)
        }else {
        gener_var_lmbds[1:9,i] <- biv_lmbds
        # Capture information about the bimodel simulation
        model_1 <- if(model[1] %in% c("rmfmkl","star")){1}else{2}
        model_2 <- if(model[2] %in% c("rmfmkl","star")){1}else{2}
        gener_var_lmbds[10:11,i]  <- c(model_1,model_2)
        }
      }
      else if(model == "rmfmkl"){
        # RMFMKL model
        gener_var_lmbds[1:5,i] <- tryCatch(c(GLDEX::fun.RMFMKL.ml(data[[i]]), 1),
                                            error = function (e){
                                              if(e$message == "non-finite value supplied by optim"){
                                                message(paste0("GLD cannot produce Lambdas with the selection of models: ",paste(model), " for variable: ",i))
                                                message(paste0("Default modgo will be used to simulate ", i))
                                                return(rep(NA, 5))
                                              }else{
                                                stop("Error in the creation of generalised lambdas")
                                              }
                                            })
      }
      else if(model == "rprs"){
        # RPRS model
        gener_var_lmbds[1:5,i] <- tryCatch(c(GLDEX::fun.RPRS.ml(data[[i]]), 2),
                                            error = function (e){
                                              if(e$message == "non-finite value supplied by optim"){
                                                message(paste0("GLD cannot produce Lambdas with the selection of models: ",paste(model), " for variable: ",i))
                                                message(paste0("Default modgo will be used to simulate ", i))
                                                return(rep(NA, 5))
                                              }else{
                                                stop("Error in the creation of generalised lambdas")
                                              }
                                            })
      }
      else if(model == "star"){
        # STAR model
        gener_var_lmbds[1:5,i] <- tryCatch(c(GLDEX::starship(data[[i]])$lambda, 1),
                                            error = function (e){
                                              if(e$message == "non-finite value supplied by optim"){
                                                message(paste0("GLD cannot produce Lambdas with the selection of models: ",paste(model), " for variable: ",i))
                                                message(paste0("Default modgo will be used to simulate ", i))
                                                return(rep(NA, 5))
                                              }else{
                                                stop("Error in the creation of generalised lambdas")
                                              }
                                            })
      }
      else if(model == "gp"){
        # Generalized Poisson distribution
        gener_var_lmbds[1:2,i] <- gp::gp.mle(data[[i]])[c("theta","lambda")]
      }
      
    }else{
      if(i %in% bin_variables){
        # Multi suggestive proportion for binary variables
        if(i %in% multi_sugg_prop){
          gener_var_lmbds[1,i] <- multi_sugg_prop[i]
          }else{
        # Calculate data set proportion
          gener_var_lmbds[1,i] <-  table(data[[i]])["1"]/sum(table(data[[i]]))
          }
      }else{
        # Default model RMFMKL
        model = "rmfmkl"
        gener_var_lmbds[1:5,i] <- tryCatch(c(GLDEX::fun.RMFMKL.ml(data[[i]]), 1),
                                                              error = function (e){
                                                                if(e$message == "non-finite value supplied by optim"){
                                                                  message(paste0("GLD cannot produce Lambdas with the selection of models: ",paste(model), " for variable: ",i))
                                                                  message(paste0("Default modgo will be used to simulate ", i))
                                                                return(rep(NA, 5))
                                                                }else{
                                                                  stop("Error in the creation of generalised lambdas")
                                                                }
                                                                })
      }
    } 
    
  }
  gener_var_lmbds <- as.data.frame(gener_var_lmbds)
  return(gener_var_lmbds)
}