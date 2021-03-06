#' Vizualize_MultAlignment (TESTING)
#'
#' Converts Multiple Alignments object to data.frame
#' @param AAS AAStringSet containing sequences to align
#' @param restrict_to FALSE=prints full alignment, "Sequence": only plot range, where sequence given in @param restrict_by is aligned. "Range", only plot range given in restrict_by
#' @param sortbyIdentity if restrict_to="Sequence", order by similarity to the Sequence named in restrict_by
#' @param DistanceMatrix if restrict_to="Sequence", optional 2x lengh of Sequene named in restrict_by matrix arranging the x position of the reads. useful for plotting alongside ab1 grahps
#' @param ROI Range of pos values. If given, marks region on plot in the defined range.
#' @importFrom reshape melt
#' @import ggplot2
#' @importFrom ggpubr theme_transparent
#' @importFrom msa msa
#' @importFrom IRanges IRanges
#' @return ggplot2 showing aligned sequences
#' @export
Vizualize_MultAlignment<-function(AAS,restrict_to=FALSE,restrict_by,sortbyIdentity=TRUE,DistanceMatrix=NULL,ROI=NULL){
  AAS<-msa(AAS)
  AAS<-DNAMultipleAlignment(AAS)
  #AAS<-methods::as(AAS,"DNAMultipleAlignment")
  AAS<-DNAMultipleAlignment_to_df(AAS)
  if(restrict_to=="Sequence"){
    AAS<-AAS[,min(which(AAS[restrict_by,]!="-")):max(which(AAS[restrict_by,]!="-"))]
  }
  if(restrict_to=="Range"){
    AAS<-AAS[,restrict_by[1]:restrict_by[2]]
  }

  if(isTRUE(sortbyIdentity) & restrict_to=="Sequence"){
    identity<-0
    for (i in 1:dim(AAS)[2]){
      identity<-identity+as.numeric(AAS[rownames(AAS)==restrict_by,i]==AAS[,i] & AAS[rownames(AAS)==restrict_by,i]!="-")
    }
    ord<-rownames(AAS[order(identity,decreasing = FALSE),])
    ident<-matrix(nrow=dim(AAS)[1],ncol=dim(AAS)[2])
    for (i in 1:dim(AAS)[2]){
      ident[,i]<-AAS[,i]==AAS[1,i]
    }
    rownames(ident)<-rownames(AAS)
  }

  AAS<-melt(AAS)
  colnames(AAS)<-c("Var1","Var2","value")


  if(isTRUE(sortbyIdentity) & restrict_to=="Sequence"){
    ident<-melt(ident)
    AAS<-cbind(AAS,ident$value)
    AAS$Var1 <- factor(AAS$Var1, levels = ord)
    colnames(AAS)[4]<-"ident"
  }

  if (!is.null(DistanceMatrix) & restrict_to=="Sequence"){
    colnames(DistanceMatrix)<-c("ID","pos")
    AAS<-merge(AAS,DistanceMatrix,by.y="ID",by.x="Var2",all.x=FALSE)
  }else{
    AAS$pos<-AAS$Var2
  }
  alignplot<-ggplot2::ggplot(AAS[AAS$Var2 & AAS$Var2,],ggplot2::aes(x=pos,y=Var1,label=value))+
    ggplot2::geom_text(size=8/ggplot2:::.pt)+
    ggplot2::scale_fill_manual(values=c("white","gray70"))+
    theme_transparent()+
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::xlab("Sequence")+
    ggplot2::theme(legend.position = "none",axis.text = ggplot2::element_text(size=ggplot2::theme_get()$text$size),axis.text.x = ggplot2::element_blank(),axis.title = ggplot2::element_text(size=8),axis.title.y = ggplot2::element_blank())

  #add ROI
  if (!is.null(ROI) & restrict_to=="Sequence"){
    if (is.null(DistanceMatrix)){
      alignplot<-addROItoGraph(alignplot,IRanges(unique(AAS$pos[AAS$Var2==ROI[1]])+0.5,unique(AAS$pos[AAS$Var2==ROI[2]])))
    } else{
      alignplot<-addROItoGraph(alignplot,IRanges(unique(AAS$pos[AAS$Var2==ROI[1]])-0.5,unique(AAS$pos[AAS$Var2==ROI[2]]+0.5)))
    }
  }
  if (!is.null(ROI) & restrict_to=="Range"){
    alignplot<-addROItoGraph(alignplot,IRanges(unique(AAS$pos[AAS$Var2==ROI[1]])-0.5,unique(AAS$pos[AAS$Var2==ROI[2]]+0.5)))
  }
  # highlight Regions
  if(restrict_to=="Sequence"){
    mark<-levels(AAS$Var1)[levels(AAS$Var1)!=restrict_by]
    for (i in mark){
      min(AAS$pos[AAS$Var1==i & AAS$value!="-"])
      alignplot<-alignplot+ggplot2::annotate("rect",
                                    xmin = min(AAS$pos[AAS$Var1==i & AAS$value!="-"])-0.5,
                                    xmax= max(AAS$pos[AAS$Var1==i & AAS$value!="-"])+0.5,
                                    ymin = which(levels(AAS$Var1)==i)-0.5, ymax = which(levels(AAS$Var1)==i)+0.5,  fill = "gray50", alpha=.3)
    }
  }



  return(alignplot)
}
