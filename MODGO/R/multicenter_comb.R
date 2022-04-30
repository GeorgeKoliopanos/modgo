#' Modgo multi-studies
#' 
#' Combines modgo objects from a multiple studies to a single one in order 
#' to calculate new correlations and visualise the data
#' 
#' 
#' @param modgo_1 a list modgo object
#' @export
#' @return A modgo object/list.
#' @author Francisco M. Ojeda, George Koliopanos
#' @keywords Multi-studies
#'
multicenter_comb <- function(modgo_1,...){
 
  argg <- c(as.list(environment()), list(...))
  
  
  nrep <- length(argg[[names(argg)[1]]][["SimulatedData"]])
  SimulatedData <- vector(mode = "list", length = nrep)
  Correlations <- vector(mode = "list", length = nrep+1)
  mean_corr <- matrix(0,
                    nrow = ncol(argg[[names(argg)[1]]][["SimulatedData"]][[1]]),
                    ncol = ncol(argg[[names(argg)[1]]][["SimulatedData"]][[1]]))
  
  for (i in c(1:nrep)){
    SimulatedData[[i]] <- argg[[names(argg)[1]]][["SimulatedData"]][[i]]
    for (j in c(2:length(names(argg)))) {
      SimulatedData[[i]] <- rbind(SimulatedData[[i]],
                                 argg[[names(argg)[j]]][["SimulatedData"]][[i]]) 
      
    }
    
    
    
    Correlations[[i]] <- cor(SimulatedData[[i]])
    mean_corr <- mean_corr + (cor(SimulatedData[[i]])/nrep)
      
  }
 
  OriginalData <- argg[[names(argg)[1]]][["OriginalData"]]
  
  for (j in c(2:length(names(argg)))) {
    OriginalData <- rbind(OriginalData,
                                argg[[names(argg)[j]]][["OriginalData"]]) 
    
  }
      Correlations[[nrep+1]] <- mean_corr
      names(Correlations) <- c(paste0("rep",seq(1:nrep)),"Mean")
      bin_variables <- argg[[names(argg)[1]]][["Binary_variables"]]
      categ_variables <- argg[[names(argg)[1]]][["Categorical_variables"]]
      combin_object <-list(SimulatedData,OriginalData,Correlations,
                           bin_variables,categ_variables)
      names(combin_object) <-c("SimulatedData","OriginalData","Correlations", 
                               "Binary_variables","Categorical_variables")
      
      return(combin_object)
      
  
  
}
