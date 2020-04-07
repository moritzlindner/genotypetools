#' OverlapAnalysis
#'
#' Make Ab1 plot
#' @param file path to Ab1 file
#' @param range limits of range to print ab1 plot
#'
#' @export
Ab1Plot<-function(file="15.4I.ab1",range, ratio = NULL){
  sangerseq <- sangerseqR::readsangerseq(file)
  # gather Ab1 traces
  traces<-as.data.frame(sangerseqR::traceMatrix(sangerseq))
  colnames(traces)<-c("A","C","G","T")
  traces$pos<-1:length(traces$A)
  traces$peak<-NA
  traces$Seq<-NA
  traces$peak[sangerseqR::peakPosMatrix(sangerseq)[,1]]<-1:length(sangerseqR::peakPosMatrix(sangerseq)[,1])
  traces$Seq[!is.na(traces$peak)]<-as.character(unlist(as.data.frame(sangerseqR::primarySeq(sangerseq))))
  `%>%`<-tidyr::`%>%`
  traces<- traces %>% tidyr::fill(peak,.direction ="up")
  traces<- traces[traces$peak %in% range[1]:range[2],]
  traces<-tidyr::gather(traces,"Base", "Amp",-c(pos,peak,Seq))
  traces$Seq[traces$Seq!=traces$Base & traces$Seq!="N"]<-NA
  ab1<-ggplot2::ggplot(data=traces,ggplot2::aes(x=pos,y=Amp,colour=Base,label=Seq))+
    ggplot2::geom_line()+ggplot2::theme_void()+
    ggplot2::geom_text(y=-100, size=8/ggplot2:::.pt)+
    ggplot2::theme(legend.position = "none")
  return(ab1)
}
