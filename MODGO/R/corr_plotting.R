#' Correlations plotting
#'
#' Plotting the correlation matrixed for Original dataset, Simulated and
#' Mean correlation matrix
#'
#' @param Modgo_obj A list object produced from modgo package
#' @param sim_dataset A number indicating the number of the simulated dataset
#' @param variables A character list listing the name of the requested variables
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
#' @importFrom ggcorrplot ggcorrplot
#' @import patchwork



corr_plots <- function(Modgo_obj,
                       sim_dataset = 1,
                       variables = colnames(Modgo_obj[["SimulatedData"]][[1]])) {
  if (!all(variables %in% colnames(Modgo_obj[["SimulatedData"]][[1]]))) {
    stop("Not all variables are in column names of data ")
    
  }
  
  
  cor1 <- as.matrix(Modgo_obj$PreSim_Sigma[variables,variables])
  cor2 <-
    cor(as.matrix(Modgo_obj[["SimulatedData"]][[sim_dataset]][, variables]))
  cor3 <- Modgo_obj[["Correlations"]][["Mean"]][variables, variables]
  cor4 <- cor1 - cor2
  
  
  original_plot <- ggcorrplot::ggcorrplot(cor1,
                                          method = "circle",
                                          title = "Original",
                                          colors = c("red","white","blue"))
  
  simulation_plot <- ggcorrplot::ggcorrplot(cor2,
                                            method = "circle",
                                            title = "Simulated",
                                            colors = c("red","white","blue"))
  
  
  mean_plot <- ggcorrplot::ggcorrplot(cor3,
                                      method = "circle",
                                      title = "Mean correlation of \nsimulations",
                                      colors = c("red","white","blue"))
  
  
  difference_plot <- ggcorrplot::ggcorrplot(cor4,
                                            method = "circle",
                                            title = "Original minus \nsimulated",
                                            colors = c("red","white","blue"))
  
  
  
  (original_plot + simulation_plot)/(mean_plot + difference_plot)
}
