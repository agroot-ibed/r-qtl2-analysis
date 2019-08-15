# QTL analysis with R/qtl2
This repository contains scripts and datasets for a QTL analysis with the [R/qtl2 package from Karl Broman](https://kbroman.org/qtl2/).  

User need to provide:
- a VCF file containing the SNP variants (previously identified through genomic analyses).
- a phenotype file containing values measured  . 

For more specific insights about the required input files, the user can consult the specific vignette: https://kbroman.org/qtl2/assets/vignettes/input_files.html 

## Datasets

Examples of dataset can be found in the `data` folder:
- **populations.snps.vcf**: a Variant Calling Format files listing the SNPs for each individual. 
- **pheno.csv**: the phenotypic values measured on each individual.
- **covar.csv**: if known covariates (such as sex) are metadata regarding the phenotypic values. 

## Usage

In R or RStudio, open and run the `qtl2_analysis.R` script to generate the required input file and perform the QTL analysis.

## Dependencies

All dependencies are installed when running the `qtl2_analysis.R` script in R or RStudio. This script makes use of:

- **VariantAnnotation**: https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html 
- **data.table**: https://cran.r-project.org/web/packages/data.table/index.html
- **R/qtl2**: https://kbroman.org/qtl2/

## Acknowledgments
We are grateful to Karl Broman for advices and support on how to use his dedicated R/qtl2 package. 

