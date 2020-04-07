#' A function to extractg the genotpe based on ab1 files
#'
#' This function retrieves ab1 files either from a link or from a directory and compares the sequence to .
#' @param dir Source directory for .ab1 files tp analyze or, if link ist not blank, directory where to download and unpack the archive containing .ab1 files to. In that case, a new subdirectory named by the current date will be created. If left blank, a file open dialog will appear.
#' @param link Optional. A url to a archive containing the ab1 files to analyze.
#' @param wtseq Sequence for wild-type
#' @param mutseq Sequence for mutant
#' @param revcomp make reverse complementary of wtseq and mutseq
#' @param cutoff minimum relative amplitude of secondary peak
#' @param separator A separator char or string by which the ab1 filenames should be split and analyzed in order to identify unique strings (e.g. mouse identifiers)
#' @return Data Frame with list of mice and respective genotypes
#' @import sangerseqR
#' @import stringr
#' @import utils
#' @import Biostrings
#' @export
sangergenotype<-function (dir = choose.dir(), link = "", wtseq = "", mutseq = "",
                          revcomp = TRUE, cutoff = 0.33, separator = "_") {
  if (revcomp == TRUE) {
    wt <- toString(Biostrings::reverseComplement(Biostrings::DNAString(wtseq)))
    mut <- toString(Biostrings::reverseComplement(Biostrings::DNAString(mutseq)))
  }
  else {
    wt <- wtseq
    mut <- mutseq
  }
  posdiff <- which(!unlist(strsplit(wt, split = "")) == unlist(strsplit(mut,
                                                                        split = "")))
  setwd(dir)
  if (!link == "") {
    dir <- paste0(dir, "\\Sequencing_", format(Sys.Date(),
                                               "%Y-%m-%d"))
    dir.create(dir)
    setwd(dir)
    utils::download.file(link, paste0(format(Sys.Date(), "%Y-%m-%d"),
                               ".", tools::file_ext(link)))
    unzip::unzip(paste0(format(Sys.Date(), "%Y-%m-%d"), ".", tools::file_ext(link)),
          junkpaths = TRUE)
  }
  files <- list.files(path = dir, pattern = "*.ab1", full.names = FALSE,
                      recursive = FALSE)
  files <- as.data.frame(files)
  files$genotype <- ""
  files$files <- as.character(files$files)
  tmp <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(files$files))
  tmp <- as.data.frame(strsplit(tmp, split = separator))
  files <- cbind(files, t(tmp[lapply(apply(tmp, 1, unique),
                                     length) > 1, ]))
  print("Extracting genotypes")
  pb <- utils::txtProgressBar(min = 0, max = length(files$files),
                       style = 3)
  for (i in 1:length(files$files)) {
    currab1 <- sangerseqR::makeBaseCalls(readsangerseq(files$files[i]),
                                         cutoff)
    allels <- c(sangerseqR::primarySeq(currab1, string = TRUE),
                sangerseqR::secondarySeq(currab1, string = TRUE))
    if (all(str_detect(allels, wt))) {
      files$genotype[i] <- "WT"
    }
    else {
      if (all(str_detect(allels, mut))) {
        files$genotype[i] <- "Hom"
      }
      else {
        if (sum(stringr::str_detect(allels, wt), stringr::str_detect(allels,
                                                                     mut)) == 2) {
          files$genotype[i] <- "Het"
        }
        else {
          files$genotype[i] <- "failed"
        }
      }
    }
    utils::setTxtProgressBar(pb, i)
  }
  rownames(files) <- 1:nrow(files)
  return(files)
}
