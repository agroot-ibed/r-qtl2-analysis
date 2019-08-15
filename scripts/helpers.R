# helpers functions
library("tidyverse")
library('VariantAnnotation')


vcf2genotypes <- function(vcfFile){
  vcf = readVcf(vcfFile)
  snp.matrix <- genotypeToSnpMatrix(vcf, uncertain = FALSE)
  snp.matrix = snp.matrix$genotypes # extract genotypes
  snp.matrix = as(snp.matrix, "character") # encodes A/A, A/B, B/B
  snp.matrix.transposed = t(snp.matrix) # to comply with the R/qtl2 genotype input format
  genotypeFile = snp.matrix.transposed 
  return(genotypeFile)
}

extract_physical_map <- function(genotypeFile){
  # from the vcf2genotypes result object, extracts a physical map file
  markers = row.names(genotypeFile)
  markers.df = data.frame(marker=markers,markers = markers)
  markers.df.parsed = markers.df %>% 
    separate(markers,into = c("chr","pos"),sep=":") 
  return(markers.df.parsed)
}
  