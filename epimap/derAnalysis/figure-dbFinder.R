## Adapted from /home/epi/ajaffe/Lieber/Projects/derfinderPaper/figure1_stem.R

## Usage:
# qrsh -l mem_free=80G,h_vmem=90G
# module load R/3.3
# mkdir -p logs
# Rscript figure-dbFinder.R > logs/figure-dbFinder_log.txt 2>&1


library('derfinder')
library('derfinderHelper')
library('derfinderPlot')
library('GenomicRanges')
library('TxDb.Hsapiens.UCSC.hg19.knownGene')
library('RColorBrewer')
library('scales')
library('GenomeInfoDb')
library("GenomicFeatures")
library('bumphunter')

## Define paths
#mainPath <- '/dcl01/lieber/ajaffe/derRuns/derSupplement/'
mainPath <- '/dcl01/lieber/ajaffe/derRuns/derChIP/epimap/'
covPath <- file.path(mainPath, 'CoverageInfo/')
resPath <- file.path(mainPath, 'derAnalysis/run1-v1.5.38-H3K4me3')
plotdir <- file.path(mainPath, 'derAnalysis', 'plots')

## Load data
message(paste(Sys.time(), 'loading results'))
load(file.path(resPath, 'colsubset.Rdata'))
load(file.path(resPath, 'groupInfo.Rdata'))
load(file.path(resPath, 'fullRegions.Rdata'))
load(file.path(resPath, 'models.Rdata'))

## Define group labels
groupSimple <- groupInfo
ageCut <- 52
levels(groupSimple) <- gsub(paste0('\\[23,', ageCut, '\\)'), paste0(ageCut, '-'), levels(groupSimple))
levels(groupSimple) <- gsub(paste0('\\[', ageCut, ',65\\]'), paste0(ageCut, '+'), levels(groupSimple))
levels(groupSimple) <- gsub('_', ':', levels(groupSimple))

## Phenotype information
message(paste(Sys.time(), 'loading phenotype information'))
load('/dcl01/lieber/ajaffe/psychENCODE_Data/EpiMap/annotated_phenotype_EpiMap_ChIPseq.rda')
pd <- pd[colsubset, ]
stopifnot(all(pd$HistoneMark == 'H3K4me3'))

## Define files
files <- pd$bamFile
names(files) <- pd$Sample_ID

## Some options
scalefac <- 32
pad <- 600

## Get best cluster of ders
ders <- fullRegions
cluster <- data.frame(area = ders$area,
    clusterChr = paste0(as.integer(ders$cluster),
    chr = as.character(seqnames(ders))))
regionClustAreas <- tapply(cluster$area, cluster$clusterChr, sum)
bestArea <- sapply(names(head(sort(regionClustAreas, decreasing=TRUE),
    20)), function(y) { which(cluster$clusterChr == y)[[1]]})
    
## Selected DER
s <- bestArea[8]
s.cluster <- range(ders[which(cluster$clusterChr == names(s))])
selected <- resize(s.cluster, width(s.cluster) + 2 * pad, fix = 'center')

## Load coverage and annotation
chr <- as.character(seqnames(selected))
chr
cov <- loadCoverage(files = files, which = selected, chr = chr)
load(file.path(resPath, chr, 'annotation.Rdata'))

## Bases
pos <- start(selected):end(selected)

## Log2 transform coverage
cov.log <- cov$coverage[pos, ]
for(i in seq_len(ncol(cov.log))) {
    cov.log[[i]] <- log2(cov.log[[i]] + scalefac)
}

## Calculate F-stats
fstats <- fstats.apply(data = cov.log, mod = models$mod, mod0 = models$mod0, scalefac = scalefac)
fstats.num <- as.numeric(fstats)
summary(fstats)


## Panel 1
bpInd <- start(ders[s]) + 100:103
ind <- start(ders[s]) - start(selected) + 100:103
covDat <- as.data.frame(cov$coverage[pos, ])
covDat.log <- as.data.frame(cov.log)
ylim <- log2(c(0, ceiling(max(covDat[ind + 1, ]))) + scalefac)
ylim[2] <- round(ylim[2] + 0.5)
y.axis <- c(0, 0.5, 2^(0:9))

group.pl <- brewer.pal(12, 'Paired')[5:12]

for(i in seq(along=ind)) {
	pdf(file.path(plotdir, paste0("cov_part", i, ".pdf")), h=6, w=4.75)
	palette(group.pl)

	y <- as.numeric(log2(covDat[ind[i] + 1,] + scalefac))
	boxplot(y ~ groupSimple, outline=FALSE, yaxt="n", main="", ylim = ylim,
        xaxt="n")
	mtext(paste0(chr, ":", bpInd[i]), line=0.5, cex=1.6)
	axis(2, at = log2(y.axis + scalefac), labels = y.axis, cex.axis = 1.3)
	points(y ~ jitter(as.numeric(groupSimple), amount=0.15),
		pch = 21, bg = as.numeric(groupSimple), cex=0.8)
    legend("topright", paste0("F=", round(fstats.num[ind[i] + 1], 2)), cex=1.3)
    if(i <= 4) {
        par(xpd = TRUE)
        legend("bottom", levels(groupSimple)[c(2 * i - 1, 2 * i)], col = group.pl[c(2 * i - 1, 2 * i)], cex=1.2, inset = -0.13, ncol = 2, bty = 'n', pch = 16)
    }	
	dev.off()
}

## panel 2
pdf(file.path(plotdir, 'fstat_panel.pdf'), h= 6,w=14)
plot(fstats.num ~ pos, type="l", xlab=chr, ylab="", cex.axis=1.4, cex.lab=1.8)
cutoff <- 4.13844212471188
abline(h=cutoff, lty=2)

pl <- rep(brewer.pal(5, 'Greys')[5], 2)
palette(pl)
sl = slice(fstats.num, lower = cutoff)
for(i in seq(along=sl)) {
	Ind = start(sl)[i]:end(sl)[i]
	polygon(x = c(pos[Ind], rev(pos[Ind])),
		y = c(fstats.num[Ind], rep(cutoff, length(Ind))),
		col = i, density = 60)
}
abline(v=range(bpInd), col="red")
dev.off()

## Smooth f-stats
fstats.smooth <- locfitByCluster(fstats.num, pos, cluster = rep(1, length(pos)), minNum = 100, bpSpan = 300, minInSpan = 100)
all(fstats.smooth$smoothed)
fstats.smooth <- fstats.smooth$fitted

## Which DERs are significant
f_ders <- ders[countOverlaps(ders, selected) > 0]
f_ders <- f_ders[order(start(f_ders), decreasing = FALSE)]

## panel 3
pdf(file.path(plotdir, 'fstat_smooth_panel.pdf'), h= 6,w=14)
plot(fstats.smooth ~ pos, type="l", xlab=chr, ylab="", cex.axis=1.4, cex.lab=1.8)
cutoff <- 4.13844212471188
abline(h=cutoff, lty=2)

pl <- rep(brewer.pal(3, 'Dark2')[1], 2)
palette(pl)
sl = slice(fstats.smooth, lower = cutoff)
for(i in seq(along=sl)) {
	Ind = start(sl)[i]:end(sl)[i]
	polygon(x = c(pos[Ind], rev(pos[Ind])),
		y = c(fstats.smooth[Ind], rep(cutoff, length(Ind))),
		col = i, density = 60)
}
abline(v=range(bpInd), col="red")
dev.off()


## coverage panel
y.axis.sample <- c(y.axis, 2^c(10:12))
pdf(file.path(plotdir, 'fullCov_panel.pdf'), h= 6,w=14)

sample.pl <- mapply(function(col, n) {
    ## Need 1/3 of lines for full saturation
    alpha(col, 1)
}, group.pl, table(groupSimple))

palette(sample.pl)
matplot(pos, covDat.log, yaxt="n",
	col=as.numeric(groupSimple), lty=1, type="l",
	xlab=chr, ylab="", cex.axis=1.4, cex.lab=1.8)
axis(2, at = log2(y.axis.sample + scalefac), labels = y.axis.sample, cex.axis = 1.3)
#m = max(covDat.log)
m <- log2(512 + scalefac)
for(i in seq(along=sl)) {
	Ind = start(sl)[i]:end(sl)[i]
	rect(xleft=min(pos[Ind]), xright = max(pos[Ind]),
		ybot = log2(scalefac), ytop =m, col=brewer.pal(5, 'Greys')[3], density=10)
}
palette(group.pl)
legend("topright", levels(groupSimple), col=seq_len(length(levels(groupSimple))), cex=1.4,pch=15, ncol = 4, bty = 'n')
dev.off()


## annotate
load("/home/epi/ajaffe/GenomicStates/GenomicState.Hsapiens.ensembl.GRCh37.p12.rda")
ensemblAnno <- annotateRegions(selected,
    GenomicState.Hsapiens.ensembl.GRCh37.p12$fullGenome)
ensemblCount <- ensemblAnno$countTable

### gene plot
a = as.data.frame(ensemblAnno$annotationList)
Strand = ifelse(a$strand == "+", 1, ifelse(a$strand=="-", -1, 0))
Col = ifelse(a$theRegion == "exon", "blue", ifelse(a$theRegion == "intron", "lightblue","white"))
Lwd = ifelse(a$theRegion == "exon", 1, ifelse(a$theRegion == "intron",0.5,0))

pdf(file.path(plotdir, 'gene_anno.pdf'), h=3,w=14)
plot(0,0, type="n", xlim=range(pos),ylim=c(-1.5,1.5),yaxt="n",ylab="",
	xlab=paste("Chromosome", mapSeqlevels(chr, 'NCBI')),	cex.axis = 1.5, cex.lab =1.8)
axis(2,c(-1,1),c("-","+"),tick=FALSE,las=1,cex.axis = 3)
abline(h=0,lty=3)
for (j in seq_len(nrow(a))) {
	polygon(c(a$start[j], a$end[j], a$end[j], a$start[j]),
	  Strand[j]/2 + c(-0.3, -0.3, 0.3, 0.3) * Lwd[j],
	  col = Col[j])
}
e <- a[a$theRegion == "exon", ]
s2 <- Strand[a$theRegion == "exon"]
g = unlist(e$symbol)
g[is.na(g)] = ""
if (length(g) > 0) {
	text(x = e$start + e$width/2, y = s2 * 0.8, g,
	  font = 2, pos = s2 + 2, cex = rep(c(1.2, 0.5, 1.2, 0.01), c(3, 3, 1, 1)))
}
dev.off()

#### extra tx info

txdb <- loadDb("/home/epi/ajaffe/Lieber/Projects/RNAseq/Ribozero_Compare/TxDb.Hsapiens.BioMart.ensembl.GRCh37.p12/inst/extdata/TxDb.Hsapiens.BioMart.ensembl.GRCh37.p12.sqlite")
txdb <- keepSeqlevels(txdb, mapSeqlevels(chr, 'NCBI'))
seqlevelsStyle(txdb) <- 'UCSC'
tx=exonsBy(txdb)
eList = tx[subjectHits(findOverlaps(selected, tx) )]

pdf(file.path(plotdir, 'exons_anno.pdf'), h=2.5,w=14)
plot(0,0, type="n", xlim=range(pos),ylim=c(0.5,length(eList)+0.5),
	yaxt="n",ylab="", xlab=paste("Chromosome", mapSeqlevels(chr, 'NCBI')), cex.axis = 1.5, cex.lab =1.8)
for(i in seq(along=eList)) {
	a = as.data.frame(eList[[i]])
	for (j in seq_len(nrow(a))) {
		polygon(c(a$start[j], a$end[j], a$end[j], a$start[j]), 
			c(i-0.25, i-0.25, i+0.25, i+0.25), col="blue")
	}
	
	int = gaps(eList[[i]])
	int = int[seqnames(int) == unique(seqnames(eList[[i]]))]
    int <- int[ end(int) < seqlengths(int) & start(int) > 1]
	end(int) = end(int)+1
	int = as.data.frame(int[start(int) != 1])
	
	for (j in seq_len(nrow(int))) {
		polygon(c(int$start[j], int$end[j], int$end[j], int$start[j]), 
			c(i-0.15, i-0.15, i+0.15, i+0.15), col="lightblue")
	}
}
dev.off()

## Reproducibility info
library('devtools')
options(width = 120)
session_info()
Sys.time()
proc.time()
