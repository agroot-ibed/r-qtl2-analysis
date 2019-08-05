vcf <- readVcf("populations.snps.vcf")
matr <- genotypeToSnpMatrix(vcf, uncertain = FALSE)
t <- t(as(matr$genotypes, "character"))
write.table(t, "SNP_matrix_m.csv", sep="\t")

