#' Ab1Plot (TESTING)
#'
#' Make Ab1 plot
#' @param file path to Ab1 file
#' @param range limits of range to print ab1 plot
#' @param ROIs Object of type GRranges with ROI(s)
#' @param ROIfill ROI fill colour
#' @importFrom sangerseqR readsangerseq traceMatrix peakPosMatrix primarySeq
#' @importFrom grid convertHeight
#' @import ggplot2 tidyr
#' @export
Ab1Plot<-function(file="15.4I.ab1",range,ROI=NULL,ROIfill="blue",size=0.1){
  sangerseq <- readsangerseq(file)
  # gather Ab1 traces
  traces<-as.data.frame(traceMatrix(sangerseq))
  colnames(traces)<-c("A","C","G","T")
  traces$pos<-1:length(traces$A)
  traces$peak<-NA
  traces$Seq<-NA
  traces$peak[peakPosMatrix(sangerseq)[,1]]<-1:length(peakPosMatrix(sangerseq)[,1])
  traces$Seq[!is.na(traces$peak)]<-as.character(unlist(as.data.frame(primarySeq(sangerseq))))
  `%>%`<-tidyr::`%>%`
  traces<- traces %>% fill(peak,.direction ="up")
  #traces<- traces[traces$peak %in% range[1]:range[2],]

  traces<-gather(traces,"Base", "Amp",-c(pos,peak,Seq))
  
  textshift<-convertHeight(ggplot2::unit(ggplot2::theme_get()$text$size*12,"pt"),"native",valueOnly=T)
  
  traces$Seq[traces$Seq!=traces$Base & traces$Seq!="N"]<-NA
  ab1<-ggplot2::ggplot(data=traces,ggplot2::aes(x=pos,y=Amp,colour=Base,label=Seq))+
    ggplot2::geom_line(size=size)+ggplot2::theme_void()+
    ggplot2::geom_text(y=textshift, size=8/ggplot2:::.pt)+
    ggplot2::theme(legend.position = "none")+
    ggplot2::coord_cartesian(xlim=range(traces$pos[traces$peak %in% range]),expand=F,ylim=c(textshift*1.5,max(traces$Amp*1.1)))
  
  if(!is.null(ROI)){
    # Translate Nucleotide position of ROI into grpah position and broaden to fully include first and last nucleotide
    narrow<-range(traces$pos[traces$peak %in% (ROI+range[1]) ])
    ROI[1]<-ROI[1]-1
    ROI[2]<-ROI[2]-1
    wide<-range(traces$pos[traces$peak %in% (ROI+range[1]) ])
    ROI[1]<-(narrow[1]+wide[1])/2
    ROI[2]<-(narrow[2]+wide[2])/2
    
    ab1<-ab1+ggplot2::annotate("rect",
                          xmin =ROI[1], 
                          xmax =ROI[2], 
                          ymin = 0, ymax = Inf,  fill = ROIfill, alpha=.3)
  }
  return(ab1)    

}