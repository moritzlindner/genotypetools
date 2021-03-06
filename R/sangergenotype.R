#' A function to interpret genotype information from ab1 files
#'
#' This function retrieves ab1 files either from a link or from a directory and compares the sequence to a reference sequence. it interprets signal amplitudes in abi traces and can thereby fidderentiate between wt, hom and het.
#' 
#' @param dir Source directory for .ab1 files to analyze or, if \var{link} is not blank, directory where to download and unpack the archive containing .ab1 files to. In that case, a new subdirectory named by the current date will be created. If left blank, a file open dialog will appear.
#' @param link Optional. A url to an archive containing the ab1 files to analyze.
#' @param wtseq Sequence snipplet caracteristic for wild-type
#' @param mutseq Sequence snipplet caracteristic for mutant
#' @param revcomp make reverse complementary of wtseq and mutseq
#' @param cutoff minimum relative amplitude of secondary peak. This is the threshold for interpreting a double peak as heterozygous.
#' @param pattern A regular expression that matches a unique pattern in the file name (like the mouse name). Could be e.g. \code{"K[0-9]{4}+"} or \code{"[0-9]{1,2}_[0-9]{1,2}[a-z]{1}+"}.
#' @return Data Frame with list of mice and respective genotypes
#' @examples
#' sangergenotype(dir="D:\\", link="http://www.x.zip",wtseq="ACTGAAAA",mutseq="ACCGAAAA", revcomp = TRUE, cutoff = 0.2)
#' @importFrom Biostrings reverseComplement DNAString
#' @importFrom sangerseqR primarySeq secondarySeq makeBaseCalls
#' @importFrom stringr str_detect
#' @export
sangergenotype<-function (dir = choose.dir(), link = "", wtseq = "", mutseq = "",
                          revcomp = TRUE, cutoff = 0.33, pattern = "K[0-9]{4}+") {
  if (revcomp == TRUE) {
    wt <- toString(Biostrings::reverseComplement(Biostrings::DNAString(wtseq)))
    mut <- toString(Biostrings::reverseComplement(Biostrings::DNAString(mutseq)))
  } else {
    wt <- wtseq
    mut <- mutseq
  }
  # posdiff <- which(!unlist(strsplit(wt, split = "")) == unlist(strsplit(mut,
  #                                                                       split = "")))
  setwd(dir)
  if (!link == "") {
    dir <- paste0(dir, "\\Sequencing_", format(Sys.Date(),
                                               "%Y-%m-%d"))
    dir.create(dir)
    setwd(dir)
    download.file(link, paste0(format(Sys.Date(), "%Y-%m-%d"),
                               ".", tools::file_ext(link)))
    unzip(paste0(format(Sys.Date(), "%Y-%m-%d"), ".", tools::file_ext(link)),
          junkpaths = TRUE)
  }
  files <- list.files(path = dir, pattern = "*.ab1", full.names = FALSE,
                      recursive = FALSE)
  files <- as.data.frame(files)
  files$genotype <- ""
  files$files <- as.character(files$files)
  files$mousename<-stringr::str_extract(files$files,pattern) # extract unique identifier from file name
  print("Extracting genotypes")
  pb <- txtProgressBar(min = 0, max = length(files$files),
                       style = 3)
  for (i in 1:length(files$files)) {
    currab1 <- makeBaseCalls(readsangerseq(files$files[i]),
                                         cutoff)
    allels <- c(primarySeq(currab1, string = TRUE),
                secondarySeq(currab1, string = TRUE))
    if (all(str_detect(allels, wt))) {
      files$genotype[i] <- "WT"
    }
    else {
      if (all(str_detect(allels, mut))) {
        files$genotype[i] <- "Hom"
      }
      else {
        if (sum(str_detect(allels, wt), str_detect(allels,mut)) == 2) {
          files$genotype[i] <- "Het"
        }
        else {
          files$genotype[i] <- "failed"
        }
      }
    }
    setTxtProgressBar(pb, i)
  }
  rownames(files) <- 1:nrow(files)
  return(files)
}
