#####
# An R script to convert VCF file into a genotype data file (with AA, AB and A/B alleles for a diploid species).
# Authors: Marc Galland
# Usage: Rscript --vanilla vcf2genotypes.R --vcf [my VCF file] --out [my genotype file ready for R/qtl2]
######


###########
# Libraries
###########

if ("qtl2" %in% installed.packages()){
  suppressMessages(library("qtl2")) # loads the library
} else {
  install.packages("qtl2", repos="http://rqtl.org/qtl2cran")
  suppressMessages(library("qtl2"))
}

if ("qtl2convert" %in% installed.packages()){
  suppressMessages(library("qtl2convert")) # loads the library
} else {
  install.packages("qtl2convert", repos="http://rqtl.org/qtl2cran")
  suppressMessages(library("qtl2"))
}

if ("data.table" %in% installed.packages()){
  suppressMessages(library("data.table")) # loads the library
} else {
  install.packages("data.table")
  suppressMessages(library("data.table"))
}

if ("VariantAnnotation" %in% installed.packages() && "snpStats" %in% installed.packages()){
  suppressMessages(library("VariantAnnotation"))
  suppressMessages(library("snpStats"))
} else {
  BiocManager::install("VariantAnnotation",update = FALSE)
  BiocManager::install('snpStats',update = FALSE)
  suppressMessages(library("VariantAnnotation"))
}

# custom functions
source("scripts/helpers.R")

####################
# Reads the VCF file 
####################
  
vcf = data.table::fread(input = "data/populations.snps.vcf",
                        sep = "\t",
                        skip = "#CHROM") # starts reading the file at this position (skips the rest)


###############################
# Extracts genotype information
###############################

# extracts genotype information from VCF file (using a custom function built around Variant Annotation package)
genotypes = vcf2genotypes("data/populations.snps.vcf")
t.genotypes = t(genotypes)
genotypes.df = as.data.frame(t.genotypes)
id=row.names(genotypes)
genotypes.final = cbind.data.frame(id,genotypes)
row.names(genotypes.final)=NULL
write.csv(x = genotypes.final,file = "data/genotypes.csv",quote = F,row.names = F)

###################################
# extracts the physical marker info
###################################
phys_markers = vcf[,1:3]
phys_markers = phys_markers[,c("ID","#CHROM","POS")]
colnames(phys_markers)=c("marker","chr","pos")
row.names(phys_markers) <- NULL
write.csv(x = phys_markers,file = "data/physical_map.csv",quote = F,row.names = F)

###################################
# File creation for R/qtl2 analysis
###################################

# create the R/qtl2 control file
write_control_file(output_file = "control.yaml",
                   crosstype = "bc",
                   geno_file = "data/genotypes.csv",
                   pmap_file = "data/physical_map.csv",
                   pheno_file = "data/pheno.csv",
                   covar_file = "data/covar.csv",
                   alleles = c("A","B"),
                   sex_covar = "sex",
                   sex_codes = list(
                     "m" = "male",
                     "f" = "female"),
                   sep = ",",
                   na.strings = "NA",
                   geno_codes =list("A/A"=1, "A/B"=2, "B/B"=3),
                   geno_transposed = TRUE,   # if TRUE, markers as rows
                   pheno_transposed = FALSE, # if TRUE, phenotypes as rows
                   covar_transposed = FALSE, # if TRUE, covariates as rows
                   overwrite = TRUE)

#######################
# Read the control file
#######################
dataset = read_cross2("control.yaml")
