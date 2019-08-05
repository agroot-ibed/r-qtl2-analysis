# qtl-analysis
## dependencies

### Install package ‘VariantAnnotation’ version 1.24.5

source("https://bioconductor.org/biocLite.R")

biocLite("VariantAnnotation")

### Install package ‘qtl’ version 1.44

install.packages("qtl")

### install other dependencies

install.packages(ggplot2)
install.packages(reshape2)
install.packages(plyr)
install.packages(qtl)
install.packages(foreach)
install.packages(doParallel)
install.packages(data.table)
install.packages(car)
install.packages(grid)
install.packages(UsingR)

## Convert vcf to SNP marix
run convert.R in the directory containing the vcf file. 

## QTL analysis
run female_qtl.R in the directory containing the SNP matrix.

