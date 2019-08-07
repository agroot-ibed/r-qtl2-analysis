#loading the neccesary dependencies
library(ggplot2)
library(reshape2)
library(plyr)
library(qtl)
library(foreach)
library(doParallel)
library(data.table)
library(car)
library(grid)
library(UsingR)



################
# Result folder
###############
resDir = dir.create("results/fem_results/")


#Reading the data

QTL <- read.cross(format="csvr", file="SNP_matrix_f_13_filtered.csv", dec='.', sep=',',genotypes =c("A/A", "A/B"), estimate.map=FALSE)

################################################### QC ###########################################################
nind(QTL)
nphe(QTL)
nchr(QTL)
totmar(QTL)
nmar(QTL)
summary(QTL)

#plot mising values
plot(QTL, auto.layout=TRUE)
plotMissing(QTL, reorder=TRUE)
par(mfrow=c(1,2),las=1)
plot(ntyped(QTL), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(QTL, "mar"), ylab="No. typed individuals",
     main="No. genotypes by marker")

#plot duplicate indiviuduals
cg <- comparegeno(QTL)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
rug(cg[lower.tri(cg)])

#plot genotype frequency
g <- pull.geno(QTL)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq)/ colSums(gfreq))
72
par(mfrow=c(1,3), las=1)
for(i in 1:3) plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "AB", "BB")[i], ylim=c(0,1))
################################################ builing genetic map #####################################################
QTL <- est.rf(QTL)
#QTL <- markerlrt(QTL)
#QTL <- markerlrt(QTL) #BIJ DRAMA DATAS
checkAlleles(QTL)
#form linkage groups
lg <- formLinkageGroups(QTL, max.rf = 0.35, min.lod = 3) #0.20 2.5
table(lg[,2])
QTL <- formLinkageGroups(QTL, max.rf = 0.35, min.lod = 3, reorgMarkers = TRUE)
#sort linkage groups
for (i in 1:382){
  QTL <- orderMarkers(QTL, chr= i)
  #pull.map(QTL, chr = i)
}
################################################## export genetic map ############################################
library(LinkageMapView)
outfile = file.path("results/", "map_15_07.pdf")
lmv.linkage.plot(QTL,outfile,mapthese=1:10, ruler= TRUE, lg.col = "lightblue1", dupnbr = TRUE)
outfile = file.path( "results/", "map_15_07.1.pdf")
lmv.linkage.plot(QTL,outfile,mapthese=11:40, ruler= TRUE, lg.col = "lightblue1", dupnbr = TRUE)
outfile = file.path( "results/", "map_15_07.1.2.pdf")
lmv.linkage.plot(QTL,outfile,mapthese=41:60, ruler= TRUE, lg.col = "lightblue1", dupnbr = TRUE)
outfile = file.path( "results/", "map_15_07.1.3.pdf")
lmv.linkage.plot(QTL,outfile,mapthese=61:80, ruler= TRUE, lg.col = "lightblue1", dupnbr = TRUE)
outfile = file.path( "results/", "map_15_07.1.4.pdf")
lmv.linkage.plot(QTL,outfile,mapthese=81:100, ruler= TRUE, lg.col = "lightblue1", dupnbr = TRUE)
outfile = file.path( "results/", "map_15_07.1.5.pdf")
lmv.linkage.plot(QTL,outfile,mapthese=101:117, ruler= TRUE, lg.col = "lightblue1", dupnbr = TRUE)
outfile = file.path( "results/", "map_15_07.1.6.pdf")
lmv.linkage.plot(QTL,outfile,mapthese=117:193, ruler= TRUE, lg.col = "lightblue1", dupnbr = TRUE)
73
png(file="mygraphic1_15.png",width=800,height=500)
plot.map(QTL, chr=1:20,alternate.chrid = TRUE)
dev.off()
png(file="mygraphic2_15.png",width=800,height=500)
plot.map(QTL, chr=21:68,alternate.chrid = TRUE)
dev.off()
png(file="mygraphic3_15.png",width=800,height=500)
plot.map(QTL, chr=61:100,alternate.chrid = TRUE)
dev.off()
png(file="mygraphic4_15.png",width=800,height=500)
plot.map(QTL, chr=101:193,alternate.chrid = TRUE)

qtlPermute <- function(..., ncores = 2, scanfun = c('scanone', 'scantwo')) {
  scanFun <- match.fun(match.arg(scanfun, c(c('scanone', 'scantwo'))))
  registerDoParallel(cores=ncores)
  times(ncores) %dopar% scanFun(...)
}

QTL <- jittermap(QTL)
#QTL analysis
family <- calc.genoprob(QTL,step=1, error.prob = 0.001)
family.logRatio <- scanone(family, pheno.col=1, method="mr")
family.logRatio.perm <- qtlPermute(cross=family, pheno.col=1, n.perm=1000, method="mr", ncores=2, scanfun="scanone")
family.logRatio.perm[1:5]
summary(family.logRatio.perm, alpha=c(0.05, 0.1, 0.2))
summary(family.logRatio, perms=family.logRatio.perm, alpha = 0.5, pvalues = TRUE)
summary(family.logRatio.perm)
siglogRatio <- summary(family.logRatio.perm ,alpha=c(0.05, 0.1, 0.2))
# Plot LOD/linkage group
png(file="all_linkage.png",width=800,height=500)
dfplot <- as.data.frame(summary(family.logRatio))
logRatio <- ggplot(dfplot, aes(x=dfplot$chr, y=dfplot$lod)) +
  geom_bar(stat="identity") + ylim(c(0,70)) + labs(title="LOD score for each linkage group") +
  geom_hline(yintercept=siglogRatio[1], linetype="dashed", colour="gray10") +
  geom_hline(yintercept=siglogRatio[2], linetype="dashed", colour="gray20") +
  xlab("linkage group") + ylab("LOD")
logRatio
dev.off()

png(file="QTL.png",width=1000,height=500)
plot(family.logRatio, ylim = (c(0,10)), alternate.chrid = TRUE, chr=1:365, show.marker.names = FALSE)
76
abline(h = siglogRatio[1],lty=5, col = "purple")
abline(h = siglogRatio[2],lty=5, col = "green")
abline(h = siglogRatio[3],lty=5, col = "red")
dev.off()
# Fine-scale mapping all linkage groups
png(file="QTL_d11.png",width=1000,height=500)
plot(family.logRatio, ylim = (c(0,10)), alternate.chrid = TRUE, chr=11, show.marker.names = TRUE)
76
abline(h = siglogRatio[1],lty=5, col = "purple")
abline(h = siglogRatio[2],lty=5, col = "green")
abline(h = siglogRatio[3],lty=5, col = "red")
dev.off()
marker.ls<-as.list(pull.map(QTL, chr = 11))
marker.ls <- as.data.frame(unlist(marker.ls))
write.table(marker.ls, file = "markerls_11.csv", sep = ",")

# Fine-scale mapping one linkage group
png(file="QTL_135.png",width=1000,height=500)
plot(family.logRatio,chr=c(135), ylim =(c(0,10)),show.marker.names= TRUE, alternate.chrid = TRUE)
abline(h = siglogRatio[1],lty=5, col = "purple")
abline(h = siglogRatio[2],lty=5, col = "green")
abline(h = siglogRatio[3],lty=5, col = "red")
dev.off()
marker.ls<-as.list(pull.map(QTL, chr = 11))
marker.ls <- as.data.frame(unlist(marker.ls))
write.table(marker.ls, file = "markerls_11.csv", sep = ",")
####################################################################################
#effect plot
png(file="effect_d11_2.png",width=500,height=500)
effectplot(QTL,pheno.col=1, mname1 = "stop",add.legend = TRUE)
z<-plotPXG(QTL,pheno.col=1, marker= "stop")
dev.off()
png(file="effect_11393_2.png",width=500,height=500)
effectplot(QTL,pheno.col=1, mname1 = "11393:53:+")
z<-plotPXG(QTL,pheno.col=1, marker= "11393:53:+" )
dev.off()
boxplot(z[,2] ~ z[,1], names=c("AA", "AB"),notch=TRUE, col=(c("white","grey")), main="d11-desaturase", xlab="Genotype",ylab="log10(16:ald/z11-16:ald)")

plot(family.logRatio,chr=c(28), ylim =(c(0,10)),show.marker.names= TRUE, alternate.chrid = TRUE)
abline(h = siglogRatio[1],lty=5, col = "purple")
abline(h = siglogRatio[2],lty=5, col = "green")
abline(h = siglogRatio[3],lty=5, col = "red")

plot(family.logRatio,chr=c(56), ylim =(c(0,10)),show.marker.names= TRUE, alternate.chrid = TRUE)
abline(h = siglogRatio[1],lty=5, col = "purple")
abline(h = siglogRatio[2],lty=5, col = "green")
abline(h = siglogRatio[3],lty=5, col = "red")


plot(family.logRatio,chr=c(77), ylim =(c(0,10)),show.marker.names= TRUE, alternate.chrid = TRUE)
abline(h = siglogRatio[1],lty=5, col = "purple")
abline(h = siglogRatio[2],lty=5, col = "green")
abline(h = siglogRatio[3],lty=5, col = "red")

plot(family.logRatio,chr=c(11,28,56,77,206,207,257), ylim =(c(0,10)),show.marker.names= TRUE, alternate.chrid = TRUE)
abline(h = siglogRatio[1],lty=5, col = "purple")
abline(h = siglogRatio[2],lty=5, col = "green")
abline(h = siglogRatio[3],lty=5, col = "red")

marker.ls<-as.list(pull.map(QTL, chr = 28))
marker.ls <- as.data.frame(unlist(marker.ls))
write.table(marker.ls, file = "markerls_28.csv", sep = ",")

marker.ls<-as.list(pull.map(QTL, chr = 56))
marker.ls <- as.data.frame(unlist(marker.ls))
write.table(marker.ls, file = "markerls_56.csv", sep = ",")

marker.ls<-as.list(pull.map(QTL, chr = 77))
marker.ls <- as.data.frame(unlist(marker.ls))
write.table(marker.ls, file = "markerls_77.csv", sep = ",")

