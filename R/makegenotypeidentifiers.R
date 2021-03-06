#' A function to extract the genotype based on ab1 files
#'
#' This function retrieves ab1 files either from a link or from a directory and compares the sequence to .
#' @param litter_list A list of litter names
#' @param size_list Corresponding list with number of pubs
#' @param alpha Make alphabetic instead of numeric identifiers
#' @examples
#' identifiers<-makegenotypeidentifiers(c("KCAW4.4","KCAW5.2"),c(8,16),alpha=TRUE)
#' @export
makegenotypeidentifiers<-function(litter_list,size_list,alpha=TRUE){
  out<-list()
  if (length(litter_list)==length(size_list)){
    for (i in 1:length(litter_list)){
      batch<-1:size_list[i]
      if (alpha){
        batch<-sapply(strsplit(paste(batch),''), function(y) paste(LETTERS[as.numeric(y)], collapse = ''))
      }
      out<-c(out,paste0(litter_list[i],batch))
    }
  }else{
    print("Unequal length of lists")
    break
  }
  unlist(out)

}
