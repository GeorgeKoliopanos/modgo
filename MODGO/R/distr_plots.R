#' Distribution plotting
#' 
#' Plotting the correlation matrices for Original dataset, Simulated and
#' Mean correlation matrix
#' 
#' @param Modgo_obj A list object produced from modgo package
#' @return A plot.
#' @author Andreas Ziegler, Francisco Miguel Echevarria, George Koliopanos
#' 
#' @examples 
#' data("Cleveland",package="modgo")
#' test_modgo <- modgo(data = Cleveland,
#'      bin_variables = c("CAD","HighFastBloodSugar","Sex","ExInducedAngina"),
#'      categ_variables =c("Chestpaintype"))
#' 
#' distr_plots(test_modgo)
#'
#' @export
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @import wesanderson

distr_plots <- function(Modgo_obj,
                        variables=colnames(Modgo_obj[["OriginalData"]]),
                        sim_dataset=1,wespalette="Cavalcanti1") {
  
  if (!all(variables %in% colnames(Modgo_obj[["OriginalData"]]))){
    
    stop("Not all variables are in column names of data ")
    
  }
  
  plotlist <- list()
  
  for (i in variables){
   comb_data <- c(Modgo_obj[["OriginalData"]][[i]],
                       Modgo_obj[["SimulatedData"]][[sim_dataset]][[i]])
    dataset <- c(rep(0,length(Modgo_obj[["OriginalData"]][[i]])),
                     rep(1,length(Modgo_obj[["SimulatedData"]][[sim_dataset]][[i]])))
    df <- as.data.frame(cbind(comb_data,dataset))
    df$dataset <- as.factor(df$dataset)
    if(i %in% Modgo_obj[["Binary_variables"]] || 
       i %in% Modgo_obj[["Categorical_variables"]]){
    p <-ggplot2::ggplot(df,aes(y=dataset)) + 
      geom_bar(aes(fill=as.factor(comb_data))) +
      scale_y_discrete(labels= c("Original","Simulated")) +
      scale_fill_manual(values = wesanderson::wes_palette(n=length(unique(comb_data)),name=wespalette,type = "continuous")) +
      theme(panel.background = element_blank(),legend.title = element_blank()) +
      ylab(label=i) +
      xlab(label="")
    } else{
     p <-ggplot2::ggplot(df,aes(x=dataset, y=comb_data,color=dataset)) + 
        theme(legend.position = "none") +
        geom_boxplot() +
        theme_bw() +
        scale_x_discrete(labels= c("Original","Simulated")) +
        scale_color_manual(values = wesanderson::wes_palette(n=length(unique(dataset)), name=wespalette,type = "continuous")) +
        ylab(label=i) +
        xlab(label="") +
        theme(panel.background = element_blank(),legend.position = "none")
      
    }
     
        plotlist[[i]] <- p
  }
 
  p <- gridExtra::grid.arrange(grobs=plotlist,ncol=ceiling(length(variables)/4))
  
  return(invisible(p))
}
