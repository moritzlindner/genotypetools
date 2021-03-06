#' addROItoGraph (TESTING)
#'
#' Adds a ROI box from IRanges to a plot
#' @param plot ggplot to append to
#' @inheritParams Ab1Plot
#' @import ggplot2
#' @importFrom BiocGenerics as.data.frame
#' @export
addROItoGraph<-function(plot,ROIs,ROIfill="blue"){
  if (!is.null(ROIs)){ # show ROIs
    for (i in 1:length(ROIs)){
      plot<-plot+ggplot2::annotate("rect",
                          xmin =
                            as.data.frame(ROIs)[i,]$start-0.5, xmax = as.data.frame(ROIs)[i,]$width+
                            as.data.frame(ROIs)[i,]$start+0.5,
                          ymin = 0, ymax = Inf,  fill = ROIfill, alpha=.3)
    }
  }

  return(plot)
}
