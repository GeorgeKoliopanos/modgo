#' Correlations plotting
#' 
#' Plotting the correlation matrixed for Original dataset, Simulated and
#' Mean correlation matrix
#' 
#' @param Modgo_obj A list object produced from modgo package
#' @return A plot.
#' @author Francisco Miguel Echevarria, George Koliopanos
#' 
#' @examples
#' 
#' data("Cleveland",package="modgo")
#' 
#' test_modgo <- modgo(data = Cleveland,
#'      bin_variables = c("CAD","HighFastBloodSugar","Sex","ExInducedAngina"),
#'      categ_variables =c("Chestpaintype"))
#'      
#' corr_plots(test_modgo)
#' 
#' @export
#' @importFrom corrplot corrplot



corr_plots <- function(Modgo_obj,sim_dataset=1,
                       variables=colnames(Modgo_obj[["OriginalData"]]),tl.cex=0.5,tl.srt=90,...) {
  
  if (!all(variables %in% colnames(Modgo_obj[["OriginalData"]]))){
    
    stop("Not all variables are in column names of data ")
    
  }
  
  opar <- par()
  par(mfrow = c(1, 4))
  cor1 <- cor(as.matrix(Modgo_obj[["OriginalData"]][,variables]))
  cor2 <- cor(as.matrix(Modgo_obj[["SimulatedData"]][[sim_dataset]][,variables]))
  cor3 <- Modgo_obj[["Correlations"]][["Mean"]][variables,variables]
  cor4 <- cor1-cor2
  
  
  corrplot::corrplot(cor1, title = "\nOriginal",tl.cex = tl.cex,
                     tl.srt=tl.srt,...)
  corrplot::corrplot(cor2, title = "\nSimulated",tl.cex = tl.cex,
                     tl.srt=tl.srt,...)
  corrplot::corrplot(cor3, title = "\nMean correlation of simulations",
                     tl.cex = tl.cex,tl.srt=tl.srt,...)
  corrplot::corrplot(cor4, title = "\nOriginal minus Simulated",tl.cex = tl.cex,
                     tl.srt=tl.srt,...)
  par(mfrow = opar[["mfrow"]])
  
  
}
  