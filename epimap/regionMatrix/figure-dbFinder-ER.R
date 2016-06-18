## Adapted from https://github.com/leekgroup/derSupplement/blob/gh-pages/figure-expressed-regions/figure-expressed-regions.R

## Usage:
# qrsh
# module load R/3.3
# mkdir -p logs
# Rscript figure-dbFinder-ER.R > logs/figure-dbFinder-ER_log.txt 2>&1

library('derfinder')
library('derfinderHelper')
library('GenomicRanges')
library('TxDb.Hsapiens.UCSC.hg19.knownGene')
library('RColorBrewer')
library('scales')
library('GenomeInfoDb')
library("GenomicFeatures")
library('bumphunter')

## Define paths
mainPath <- '/dcl01/lieber/ajaffe/derRuns/derChIP/epimap/'
covPath <- file.path(mainPath, 'CoverageInfo/')
resPath <- file.path(mainPath, 'derAnalysis/run1-v1.5.38-H3K4me3')
plotdir <- file.path(mainPath, 'regionMatrix', 'figure-ER')
dir.create(plotdir, showWarnings = FALSE)


## Load data
load(file.path(mainPath, 'regionMatrix', 'regionMat-H3K4me3-cut10-chr12.Rdata'))
load(file.path(resPath, 'groupInfo.Rdata'))
load(file.path(resPath, 'models.Rdata'))
load(file.path(resPath, 'colsubset.Rdata'))

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
pad <- 300
scalefac <- 1

## Selected region
selected <- range(regionMat$chr12$regions[subjectHits(findOverlaps(GRanges('chr12', IRanges(11652500, 11654500), '*'), regionMat$chr12$regions))])
selected <- resize(selected, width(selected) + 2 * pad, fix = 'center')

## Load coverage
chr <- as.character(seqnames(selected))
chr
cov <- loadCoverage(files = files, which = selected, chr = chr, protectWhich = 3e4, totalMapped = pd$totalMapped)


## Bases
pos <- start(selected):end(selected)

## Log2 transform coverage
cov.log <- cov$coverage[pos, ]
for(i in seq_len(ncol(cov.log))) {
    cov.log[[i]] <- log2(cov.log[[i]] + scalefac)
}

## Misc
covDat <- as.data.frame(cov$coverage[pos, ])
covDat.log <- as.data.frame(cov.log)


## Calculate overall mean
mean.ov <- log2(rowMeans(covDat) + scalefac)
y.axis <- c(0, 2^(1:9))

## Mean panel
pdf(file.path(plotdir, "mean_panel.pdf"), h= 6,w=14)
plot(mean.ov ~ pos, type="l", xlab=chr, ylab="", cex.axis=1.4, cex.lab=1.8, yaxt="n", ylim = log2(c(4, max(y.axis) + scalefac)))
axis(2, at = log2(y.axis + scalefac), labels = y.axis, cex.axis = 1.2)
mean.cutoff <- log2(10 + scalefac + 1)
abline(h= mean.cutoff, lty=2)

mean.sl <- slice(mean.ov, lower = mean.cutoff)
pl <- rep(brewer.pal(5, 'Greys')[5], 2)
palette(pl)
for(i in seq(along = mean.sl)) {
	Ind = start(mean.sl)[i]:end(mean.sl)[i]
	polygon(x = c(pos[Ind], rev(pos[Ind])),
		y = c(mean.ov[Ind], rep(mean.cutoff, length(Ind))),
		col = i, density =60)
}
dev.off()

## Smooth mean panel
mean.ov.smooth <- log2(runmedByCluster(rowMeans(covDat), k = 299, cluster = rep(1, nrow(covDat)), x = pos)$fitted[, 1] + scalefac)
y.axis <- c(0, 2^(1:8))

pdf(file.path(plotdir, "mean_smooth_panel.pdf"), h= 6,w=14)
plot(mean.ov.smooth ~ pos, type="l", xlab=chr, ylab="", cex.axis=1.4, cex.lab=1.8, yaxt="n", ylim = log2(c(4, max(y.axis) + scalefac)))
axis(2, at = log2(y.axis + scalefac), labels = y.axis, cex.axis = 1.2)
abline(h= mean.cutoff, lty=2)

mean.sl <- slice(mean.ov.smooth, lower = mean.cutoff)
pl <- rep(brewer.pal(3, 'Dark2')[1], 2)
palette(pl)
for(i in seq(along = mean.sl)) {
	Ind = start(mean.sl)[i]:end(mean.sl)[i]
	polygon(x = c(pos[Ind], rev(pos[Ind])),
		y = c(mean.ov.smooth[Ind], rep(mean.cutoff, length(Ind))),
		col = i, density =60)
}
dev.off()




## coverage panel
y.axis.sample <- c(0, 2^(1:12))
group.pl <- brewer.pal(12, 'Paired')[5:12]


pdf(file.path(plotdir, "fullCov_panel.pdf"), h= 6,w=14)

sample.pl <- mapply(function(col, n) {
    alpha(col, 1)
}, group.pl, table(groupSimple))

palette(sample.pl)
matplot(pos, covDat.log, yaxt="n",
	col=as.numeric(groupSimple), lty=1, type="l",
	xlab=chr, ylab="", cex.axis=1.4, cex.lab=1.8)
axis(2, at = log2(y.axis.sample + scalefac), labels = y.axis.sample, cex.axis = 1.3)
#m = max(covDat.log)
m <- log2(512 + scalefac)
for(i in seq(along=mean.sl)) {
	Ind = start(mean.sl)[i]:end(mean.sl)[i]
	rect(xleft=min(pos[Ind]), xright = max(pos[Ind]),
		ybot = log2(scalefac), ytop = m, col=brewer.pal(5, 'Greys')[3], density=10)
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

pdf(file.path(plotdir, "gene_anno.pdf"), h=3,w=14)
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
	  font = 2, pos = s2 + 2,
      cex = c(1.2, 0.01, 0.5, 0.5, 0.5, 0.01, 1.2, 1.2, 1.2, 0.01, 0.01))
}
dev.off()

#### extra tx info

txdb <- loadDb("/home/epi/ajaffe/Lieber/Projects/RNAseq/Ribozero_Compare/TxDb.Hsapiens.BioMart.ensembl.GRCh37.p12/inst/extdata/TxDb.Hsapiens.BioMart.ensembl.GRCh37.p12.sqlite")
txdb <- keepSeqlevels(txdb, mapSeqlevels(chr, 'NCBI'))
seqlevelsStyle(txdb) <- 'UCSC'
tx=exonsBy(txdb)
eList = tx[subjectHits(findOverlaps(selected, tx) )]

e.strand <- unlist(unique(strand(eList)))
e.n.neg <- sum(e.strand == '-')
e.n.pos <- sum(e.strand == '+')
ylim <- c(-1 * e.n.neg + ifelse(e.n.neg > 0, -0.5, 0.5), e.n.pos + 0.5)

pdf(file.path(plotdir, "trans_anno.pdf"), h=4.5,w=14)
plot(0,0, type="n", xlim=range(pos), ylim=ylim,
	yaxt="n",ylab="", xlab=paste("Chromosome", mapSeqlevels(chr, 'NCBI'), '(161.1 mb)'), xaxt='n', cex.lab = 1.8)
axis(1, at = c(161115000, 161120000, 161125000, 161130000), labels = c('+15k', '+20k', '+25k', '+30k'), cex.axis = 1.5)
axis(2, c(- ifelse(e.n.neg, median(seq_len(e.n.neg)), NA), ifelse(e.n.pos, median(seq_len(e.n.pos)), NA)), c(ifelse(e.n.neg, '-', NA), ifelse(e.n.pos, "+", NA)), tick=FALSE,las=1,cex.axis = 3)
abline(h=0,lty=3)
for(i in seq(along=eList)) {
	a = as.data.frame(eList[[i]])
    i.strand <- sum(e.strand[ seq_len(length(e.strand)) <= i] == e.strand[i]) * ifelse(e.strand[i] == "+", 1, -1)
	for (j in seq_len(nrow(a))) {
		polygon(c(a$start[j], a$end[j], a$end[j], a$start[j]), 
			c(i.strand - 0.25, i.strand -0.25, i.strand +0.25, i.strand +0.25), col="blue")
	}
	
	int = gaps(eList[[i]])
	int = int[seqnames(int) == unique(seqnames(eList[[i]]))]
    int <- int[ end(int) < seqlengths(int) & start(int) > 1]
	end(int) = end(int)+1
	int = as.data.frame(int[start(int) != 1])
	
    
	for (j in seq_len(nrow(int))) {
		polygon(c(int$start[j], int$end[j], int$end[j], int$start[j]), 
			c(i.strand - 0.15, i.strand -0.15, i.strand + 0.15, i.strand +0.15), col="lightblue")
	}
    
}
dev.off()


