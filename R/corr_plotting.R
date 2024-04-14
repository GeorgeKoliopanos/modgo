#' Plots correlation matrix of original and simulated data
#'
#' Produces a graphical display of the Pearson correlation matrix of the 
#' original dataset, a  single simulated  dataset and also of the average of 
#' the correlation matrices across all simulations for an object returned by 
#' \code{\link[modgo]{modgo}}.
#' 
#' @param Modgo_obj An object returned by \code{\link[modgo]{modgo}}.
#' @param sim_dataset Number indicating the simulated dataset in 
#' \code{Modgo_obj} to be used in plots. 
#' @param variables A character vector indicating the columns in the data to 
#' be used in plots.
#' @return A patchwork object created by 
#' \code{\link[patchwork:wrap_plots]{patchwork::wrap_plots}}.
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
                       variables = colnames(Modgo_obj[["simulated_data"]][[1]])) {
  if (!all(variables %in% colnames(Modgo_obj[["simulated_data"]][[1]]))) {
    stop("Not all variables are in column names of data ")
    
  }
  
  cor1 <- as.matrix(Modgo_obj$presim_sigma[variables,variables])
  cor2 <-
    cor(as.matrix(Modgo_obj[["simulated_data"]][[sim_dataset]][, variables]))
  cor3 <- Modgo_obj[["correlations"]][["Mean"]][variables, variables]
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
