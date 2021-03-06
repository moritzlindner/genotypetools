# sgtR
### A package to accelerate interpretaition of sanger sequencing based genotyping results

## Installation
```{r}
if (!requireNamespace("remotes", quietly = TRUE)){
  install.packages("remotes")
}
remotes::install_github("moritzlindner/genotypetools")
```

## Description
The main function of this package is \code{sangergenotype()}. It takes a directory or a url to a zip file, cycles through all ab1 files therein and compares the traces to provided reference sequences for wt and mutant alleles. Based on thresholding it can discriminate between wt, het, and hom. Results are returned in a data frame. There are several additional functions which can be helpful when visualizing sanger genotyping data.

## Examples

### Example 1:

```{r Example1, eval=FALSE, include=T}
sangergenotype(dir="D:\\", link="http://www.x.zip",wtseq="ACTGAAAA",mutseq="ACCGAAAA", revcomp = TRUE, cutoff = 0.2)
as.data.frame(tmp)
```

Developed by [Moritz Lindner](https://www.uni-marburg.de/en/fb20/departments/physiology/research/dominik-oliver-lab/research2/retinal-physiology-and-gene-therapy)
