#' DNAMultipleAlignment_to_df (TESTING)
#'
#' This function is designed to convert Multiple Alignments of the XStringSet class into a data.frame
#' @param DNAS_align Multiple alignment object
#' @return DF a data.frame containing the sequences. Rows represent individual sequences, columns the sequence "items" (i.e. one base per column).
#' @importFrom Biostrings DNAStringSet
#' @export
DNAMultipleAlignment_to_df<-function(DNAS_align){
  tmp<-NULL
  for ( i in 1:length(names(Biostrings::DNAStringSet(DNAS_align@unmasked)))){
    tmp$names[i]<-names(Biostrings::DNAStringSet(DNAS_align@unmasked))[i]
    tmp$seq[i]<-as.character(Biostrings::DNAStringSet(DNAS_align@unmasked[[i]]))
  }
  DF<-array(character(),dim=c(length(tmp$seq),mean(nchar(tmp$seq[1]))))
  for (i in 1:length(tmp$seq)){
    DF[i,]<-substring(tmp$seq[i], seq(1, mean(nchar(tmp$seq)), 1), seq(1, mean(nchar(tmp$seq)), 1))
  }
  rownames(DF)<-names(Biostrings::DNAStringSet(DNAS_align@unmasked))
  return(DF)
}
