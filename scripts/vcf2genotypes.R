#####
# An R script to convert VCF file into a genotype data file (with AA, AB and A/B alleles for a diploid species).
# Authors: Marc Galland
# Usage: Rscript --vanilla vcf2genotypes.R --vcf [my VCF file] --out [my genotype file ready for R/qtl2]
######


###########
# Libraries
###########

if ("qtl2" %in% installed.packages()){
  library("qtl2") # loads the library
} else {
  install.packages("qtl2", repos="http://rqtl.org/qtl2cran")
  library("qtl2")
}

if ("qtl2convert" %in% installed.packages()){
  library("qtl2convert") # loads the library
} else {
  install.packages("qtl2convert", repos="http://rqtl.org/qtl2cran")
  library("qtl2")
}

if ("data.table" %in% installed.packages()){
  library("data.table") # loads the library
} else {
  install.packages("data.table")
  library("data.table")
}

if ("VariantAnnotation" %in% installed.packages() && "snpStats" %in% installed.packages()){
  library("VariantAnnotation")
  library("snpStats")
} else {
  BiocManager::install("VariantAnnotation",update = FALSE)
  BiocManager::install('snpStats',update = FALSE)
  library("VariantAnnotation")
}

# custom functions
source("scripts/helpers.R")

####################
# Reads the VCF file 
####################

vcf = data.table::fread(input = "variants.vcf",
                        sep = "\t",
                        skip = "#CHROM") # starts reading the file at this position (skips the rest)

#################################################
# Conversion to genotypes with the R/qtl2 package
#################################################

# I need to separate the map information from the genotype

###################################
# File creation for R/qtl2 analysis
###################################

# extracts genotype information from VCF file (Variant Annotation package)
genotypes = vcf2genotypes("variants.vcf")
write.csv(x = genotypes,file = "rqtl2_inputfiles/genotypes.csv",quote = F,row.names = T)

# extracts the physical marker info
phys_markers = extract_physical_map(genotypes)
write.csv(x = phys_markers,file = "rqtl2_inputfiles/physical_map.csv",quote = F,row.names = T)

# make a fake table for phenotypic values
pheno = data.frame(
  id = colnames(genotypes),
  pheno.vals = runif(n = length(colnames(genotypes)),min = 0.2,max = 10)
  )

# create the R/qtl2 control file
write_control_file(output_file = "rqtl2_inputfiles/control.yaml",
                   crosstype = "f2",
                   geno_file = "rqtl2_inputfiles/genotypes.csv",
                   pmap_file = "rqtl2_inputfiles/physical_map.csv",
                   pheno_file = "rqtl2_inputfiles/pheno.csv",
                   alleles = c("A","B"),
                   geno_codes = c(1,2,3),
                   sep = ",",
                   na.strings = "NA",
                   geno_transposed = TRUE,
                   overwrite = TRUE)
