#' Correlations plotting
#' 
#' Plotting the correlation matrixed for Original dataset, Simulated and
#' Mean correlation matrix
#' 
#' @param Modgo_obj A list object produced from modgo package
#' @return A plot.
#' @author Francisco M. Ojeda, George Koliopanos
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
                       variables=colnames(Modgo_obj[["OriginalData"]]),
                       tl.cex=0.5,tl.srt=90,title_pos=5,title_height=0.0,
                       title_size=0.8,...) {
  
  if (!all(variables %in% colnames(Modgo_obj[["OriginalData"]]))){
    
    stop("Not all variables are in column names of data ")
    
  }
  
  opar <- par()
  par(mfrow = c(2, 2))
  cor1 <- cor(as.matrix(Modgo_obj[["OriginalData"]][,variables]))
  cor2 <- cor(as.matrix(Modgo_obj[["SimulatedData"]][[sim_dataset]][,variables]))
  cor3 <- Modgo_obj[["Correlations"]][["Mean"]][variables,variables]
  cor4 <- cor1-cor2
  
  
  corrplot::corrplot(cor1,tl.cex = tl.cex,
                     tl.srt=tl.srt,tl.col = "black",...)
  mtext("Original", at=title_pos, line=title_height, cex=title_size)
  corrplot::corrplot(cor2,tl.cex = tl.cex,
                     tl.srt=tl.srt,tl.col = "black",...)
  mtext("Simulated", at=title_pos, line=title_height, cex=title_size)
  
  corrplot::corrplot(cor3,
                     tl.cex = tl.cex,tl.srt=tl.srt,tl.col = "black",...)
  mtext("Mean correlation of \nsimulations", at=title_pos, line=title_height, cex=title_size)
  
  corrplot::corrplot(cor4,tl.cex = tl.cex,
                     tl.srt=tl.srt,tl.col = "black",...)
  mtext("Original minus \nsimulated", at=title_pos, line=title_height, cex=title_size)
  
  par(mfrow = opar[["mfrow"]])
  
  
}
  