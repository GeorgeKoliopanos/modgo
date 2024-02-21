#' Plots distribution of original and simulated data
#' 
#' Produces a graphical display of the distribution of the variables
#' of the original dataset and a single simulated dataset for an object 
#' returned by \code{\link[modgo]{modgo}}.
#' 
#' For continuous variables box-and-whisker plots are displayed, while 
#' categorical variables bar charts are produced.
#' 
#' @inheritParams corr_plots
#' @param  wespalette a name of the selected wesanderson color pallet
#' @param  text_size a number for the  size of the annotation text
#' @return a plot.
#' @author Andreas Ziegler, Francisco M. Ojeda, George Koliopanos
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
                        variables=colnames(Modgo_obj[["original_data"]]),
                        sim_dataset=1,
                        wespalette="Cavalcanti1",
                        text_size=12) {
  
  if(is.null(Modgo_obj[["original_data"]])){
    stop("Distr_plots cannot run without providing an original dataset in the main modgo function")
  }
  if (!all(variables %in% colnames(Modgo_obj[["original_data"]]))){
    
    stop("Not all variables are in column names of data ")
    
  }
  
  plotlist <- list()
  
  for (i in variables){
   comb_data <- c(Modgo_obj[["original_data"]][[i]],
                       Modgo_obj[["simulated_data"]][[sim_dataset]][[i]])
    dataset <- c(rep(1,length(Modgo_obj[["original_data"]][[i]])),
                     rep(0,length(Modgo_obj[["simulated_data"]][[sim_dataset]][[i]])))
    df <- as.data.frame(cbind(comb_data,dataset))
    df$dataset <- as.factor(df$dataset)
    if(i %in% Modgo_obj[["bin_variables"]] || 
       i %in% Modgo_obj[["categ_variables"]]){
    p <-ggplot2::ggplot(df,aes(y=dataset)) + 
      geom_bar(aes(fill=as.factor(comb_data))) +
      scale_y_discrete(labels= c("Simulated","Original")) +
      scale_fill_manual(values = wesanderson::wes_palette(n=length(unique(comb_data)),name=wespalette,type = "continuous")) +
      theme(panel.background = element_blank(),legend.title = element_blank(),
            axis.text.y = element_text(color="black", size=text_size+2, face="bold"),
            axis.text.x = element_text(color="black", size=text_size),
            axis.title.y = element_text(color="black", size=text_size+2, face="bold"),
            legend.text = element_text(color="black", size=text_size+2,face="bold")) +
      ylab(label=i) +
      xlab(label="")
    } else{
      comb_data <- c(Modgo_obj[["original_data"]][[i]],
                     Modgo_obj[["simulated_data"]][[sim_dataset]][[i]])
      dataset <- c(rep(0,length(Modgo_obj[["original_data"]][[i]])),
                   rep(1,length(Modgo_obj[["simulated_data"]][[sim_dataset]][[i]])))
      df <- as.data.frame(cbind(comb_data,dataset))
      df$dataset <- as.factor(df$dataset)
     p <-ggplot2::ggplot(df,aes(x=dataset, y=comb_data,color=dataset)) + 
        theme(legend.position = "none") +
        geom_boxplot() +
        theme_bw() +
        scale_x_discrete(labels= c("Original","Simulated")) +
        scale_color_manual(values = wesanderson::wes_palette(n=length(unique(dataset)), name=wespalette,type = "continuous")) +
        ylab(label=i) +
        xlab(label="") +
        theme(panel.background = element_blank(),legend.position = "none",
              axis.text.x = element_text(color="black", size=text_size+2, face="bold"),
              axis.text.y = element_text(color="black", size=text_size),
              axis.title.y = element_text(color="black", size=text_size+2, face="bold"))
        
      
    }
     
        plotlist[[i]] <- p
  }
 
  p <- gridExtra::grid.arrange(grobs=plotlist,ncol=ceiling(length(variables)/4))
  
  return(invisible(p))
}
