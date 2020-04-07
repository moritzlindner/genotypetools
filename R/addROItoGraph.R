#' addROItoGraph
#'
#' Adds a ROI box from IRanges to a plot
#' @param plot ggplot to append to
#' @param ROIs Object of type GRranges with ROI(s)
#' @param ROIfill ROI fill colour
#' @import ggplot2
#' @import BiocGenerics
#' @export
addROItoGraph<-function(plot,ROIs,ROIfill="blue"){
  if (!is.null(ROIs)){ # show ROIs
    for (i in 1:length(ROIs)){
      plot<-plot+ggplot2::annotate("rect",
                          xmin =
                            BiocGenerics::as.data.frame(ROIs)[i,]$start-0.5, xmax = as.data.frame(ROIs)[i,]$width+
                            BiocGenerics::as.data.frame(ROIs)[i,]$start+0.5,
                          ymin = 0, ymax = Inf,  fill = ROIfill, alpha=.3)
    }
  }

  return(plot)
}
