#' Generalized Lambda and Poisson preparation
#' 
#' Prepare the four moments matrix for GLD and GPD
#' @param data a data frame with original variables. 
#' @param gener_var data frame with transformed variables. 
#' @param gener_var_model A numeric square matrix.
#' @return A numeric matrix
#' @export
generalizedMatrix <- function(data, gener_var, gener_var_model){
  
  # Prepare gener_var four moments
  gener_var_lmbds <- matrix(nrow = 11, ncol = length(gener_var))
  colnames(gener_var_lmbds) <- gener_var
  
  for (i in gener_var){
    if (i %in% gener_var_model[,"Variables"]){
      
      model <- unlist(strsplit(
        gener_var_model[which(gener_var_model[,"Variables"] == i), 
                        "Model"],
        split = "-"))
      if(length(model) == 2){
        biv_lmbds <- tryCatch(GLDEX::fun.auto.bimodal.ml(data[[i]],
                                                                      init1.sel=model[1],
                                                                      init2.sel=model[2],
                                                                      init1=c(-0.25,1.5),
                                                                      init2=c(-0.25,1.5),
                                                                      leap1=3,leap2=3)$par
                                           ,
                                           error = function (e){
                                             paste0("GLD cannot produce Lambdas with the selection of models: ",paste(model,collapse = " - "))
                                           }
        )
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
        gener_var_lmbds[1:4,i] <- GLDEX::fun.RMFMKL.ml(data[[i]])
        gener_var_lmbds[5,i] <- 1
      }
      else if(model == "rprs"){
        gener_var_lmbds[1:4,i] <- GLDEX::fun.RPRS.ml(data[[i]])
        gener_var_lmbds[5,i] <- 2
      }
      else if(model == "star"){
        gener_var_lmbds[1:4,i] <- GLDEX::starship(data[[i]])$lambda
        gener_var_lmbds[5,i] <- 1
      }
      else if(model == "gp"){
        gener_var_lmbds[1:2,i] <- gp::gp.mle(data[[i]])[c("theta","lambda")]
      }
      
    }else{
      gener_var_lmbds[1:4,i] <- GLDEX::fun.RMFMKL.ml(data[[i]])
      gener_var_lmbds[5,i] <- 1
    } 
    
  }
  gener_var_lmbds <- as.data.frame(gener_var_lmbds)
  return(gener_var_lmbds)
}