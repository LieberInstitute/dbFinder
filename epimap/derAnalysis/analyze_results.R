## Analyze results

## Load libraries
library('getopt')

## Available at http://www.bioconductor.org/packages/release/bioc/html/derfinder.html
library('derfinder')
library('derfinderPlot')
library('GenomicRanges')
library('limma')
library('GenomicFeatures')
library('biomaRt')
library('bumphunter')
library('RColorBrewer')
library('devtools')

## Specify parameters
spec <- matrix(c(
    'histone', 'i', 1, 'character', 'For epimap, the histone mark to use. Either H3K27ac or H3K4me3.',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## For testing
if(FALSE) opt <- list(histone = 'H3K4me3')

## Dirs
rootdir <- '/dcl01/lieber/ajaffe/derRuns/derChIP/epimap/'
maindir <- file.path(rootdir, 'derAnalysis', paste0('run1-v1.5.38-',
    opt$histone))
plotdir <- file.path(rootdir, 'derAnalysis', 'plots')
dir.create(plotdir, showWarnings = FALSE)
    
## Results
load(file.path(maindir, 'groupInfo.Rdata'))
load(file.path(maindir, 'colsubset.Rdata'))
## Load regions
message(paste(Sys.time(), 'loading fullRegions.Rdata'))
load(file.path(maindir, 'fullRegions.Rdata'))

#keepIndex <- width(fullRegions) >=6
keepIndex <- fullRegions$significantQval == 'TRUE'
regions <- fullRegions[keepIndex]

## Phenotype information
message(paste(Sys.time(), 'loading phenotype information'))
load('/dcl01/lieber/ajaffe/psychENCODE_Data/EpiMap/annotated_phenotype_EpiMap_ChIPseq.rda')
pd <- pd[colsubset, ]
stopifnot(all(pd$HistoneMark == opt$histone))

## By ENSEMBL annotation
message(paste(Sys.time(), 'loading ENSEMBL genomic state'))
load("/home/epi/ajaffe/GenomicStates/GenomicState.Hsapiens.ensembl.GRCh37.p12.rda")
gs <- GenomicState.Hsapiens.ensembl.GRCh37.p12$fullGenome

message(paste(Sys.time(), 'annotating with ENSEMBL v75'))
ensemblAnno <- annotateRegions(regions, gs)
save(ensemblAnno, file = file.path(maindir, 'ensemblAnno.Rdata'))
countTable <- ensemblAnno$countTable

## annotation ####
dim(countTable)
annoClassList <- list(strictExonic = 
	which(countTable[, "exon"] > 0 & countTable[, "intron"] == 0 &
		countTable[,"intergenic"] == 0),
	strictIntronic = 
	which(countTable[, "intron"] > 0 & countTable[, "exon"] == 0 &
		countTable[, "intergenic"] == 0),
	strictIntergenic = which(countTable[, "intergenic"] > 0 &
        countTable[, "exon"] == 0 & countTable[, "intron"] == 0),
	exonIntron = which(countTable[, "exon"] > 0 & countTable[, "intron"] > 0 &
		countTable[, "intergenic"] == 0))
annoClassList$All <- seq_len(length(regions))
print('Annotation breakdown')
sapply(annoClassList, length)
100 * sapply(annoClassList, length) / nrow(countTable)

## Region width
print('Region width information')
quantile(width(regions))
sapply(annoClassList, function(ii) quantile(width(regions[ii])))

## Venn diagram: code modified from limma::vennDiagram
vennDiagram_custom <- function (object, include = "both", names = NULL, 
    mar = rep(1, 4), cex = c(1.5, 1, 0.7), lwd = 1, circle.col = NULL,
    counts.col = NULL, text.col = NULL, ...) 
{
    include <- as.character(include)
    LenInc <- min(length(include), 2)
    if (is(object, "VennCounts")) {
        include <- include[1]
        LenInc <- 1
    }
    else {
        if (LenInc > 1) 
            z2 <- vennCounts(object, include = include[2])[, 
                "Counts"]
        object <- vennCounts(object, include = include[1])
    }
    z <- object[, "Counts"]
    nsets <- ncol(object) - 1
    if (nsets > 5) 
        stop("Can't plot Venn diagram for more than 5 sets")
    VennZone <- object[, 1:nsets, drop = FALSE]
    VennZone <- apply(VennZone, 1, function(x) paste(x, sep = "", 
        collapse = ""))
    names(z) <- VennZone
    if (length(include) == 2) 
        names(z2) <- VennZone
    if (is.null(names)) 
        names <- colnames(object)[1:nsets]
    FILL.COL <- TRUE
    if (is.null(circle.col)) {
        circle.col <- par("col")
        FILL.COL <- FALSE
    }
    if (length(circle.col) < nsets) 
        circle.col <- rep(circle.col, length.out = nsets)
    if (is.null(counts.col)) 
        counts.col <- par("col")
    if (length(counts.col) < LenInc) 
        counts.col <- rep(counts.col, length.out = LenInc)
    if(is.null(text.col)) text.col <- rep('black', switch(nsets, counts.col[1], counts.col[1], 8))
    old.par <- par()$mar
    on.exit(par(mar = old.par))
    par(mar = mar)
    if (nsets <= 3) {
        plot(x = 0, y = 0, type = "n", xlim = c(-4, 4), ylim = c(-4, 
            4), xlab = "", ylab = "", axes = FALSE, ...)
            
        theta <- 2 * pi * (0:360)/360
        xcentres <- switch(nsets, 0, c(-1, 1), c(-1, 1, 0))
        ycentres <- switch(nsets, 0, c(0, 0), c(1, 1, -2)/sqrt(3))
        r <- 2
        xtext <- switch(nsets, -1.2, c(-1.2, 1.2), c(-1.2, 1.2, 
            0))
        ytext <- switch(nsets, 1.8, c(1.8, 1.8), c(3, 3, 
            -3.5))
        for (circle in 1:nsets) {
            if (!FILL.COL) 
                lines(xcentres[circle] + r * cos(theta), ycentres[circle] + 
                  r * sin(theta), lwd = lwd, col = circle.col[circle])
            if (FILL.COL) {
                RGB <- col2rgb(circle.col[circle])/255
                ALPHA <- 0.06
                RGB.ALP <- rgb(RGB[1, 1], RGB[2, 1], RGB[3, 1], 
                  alpha = ALPHA)
                polygon(xcentres[circle] + r * cos(theta), ycentres[circle] + 
                  r * sin(theta), border = circle.col[circle], 
                  lwd = lwd, col = RGB.ALP)
            }
            text(xtext[circle], ytext[circle], names[circle], 
                cex = cex * 1.3, col = circle.col[circle])
        }
        switch(nsets, rect(-3, -2.5, 3, 2.5), rect(-3, -2.5, 
            3, 2.5), rect(-3.9, -3.9, 3.9, 3.9))
        showCounts <- switch(nsets, function(counts, cex, adj, 
            col, leg) {
            text(2.3, -2.1, counts[1], cex = cex, col = col, 
                adj = adj)
            text(0, 0, counts[2], cex = cex, col = col, adj = adj)
        }, function(counts, cex, adj, col, leg) {
            text(2.3, -2.1, counts[1], cex = cex, col = col, 
                adj = adj)
            text(1.5, 0.1, counts[2], cex = cex, col = col, adj = adj)
            text(-1.5, 0.1, counts[3], cex = cex, col = col, 
                adj = adj)
            text(0, 0.1, counts[4], cex = cex, col = col, adj = adj)
        }, function(counts, cex, adj, col, leg) {
            text(3, -3, counts[1], cex = cex, col = col[1], adj = adj)
            text(0, -2.2, counts[2], cex = cex * 1.5, col = col[2], adj = adj)
            text(2, 1, counts[3], cex = cex * 1.5, col = col[3], adj = adj)
            text(1.3, -0.5, counts[4], cex = cex, col = col[4], 
                adj = adj)
            text(-2, 1, counts[5], cex = cex * 1.5, col = col[5], adj = adj)
            text(-1.3, -0.5, counts[6], cex = cex * 1.3, col = col[6], 
                adj = adj)
            text(0, 1.3, counts[7], cex = cex, col = col[7], adj = adj)
            text(0, 0, counts[8], cex = cex, col = col[8], adj = adj)
        })
        if (LenInc == 1) 
            adj <- c(0.5, 0.5)
        else adj <- c(0.5, 0)
        showCounts(counts = z, cex = cex[1], adj = adj, col = text.col, 
            leg = include[1])
        return(invisible())
    }
}

pdf(file.path(plotdir, paste0(opt$histone, '_ensemblVenn.pdf')), width = 10,
    height = 10)
venn_col <- brewer.pal(7, "Dark2")[c(1, 3, 2, 4)]
vennDiagram_custom(vennCounts(countTable > 0), 
    main = 'dbPeaks overlap with Ensembl v75 features', cex.main = 2,
    circle.col = venn_col[1:3], lwd = 1.5, cex = 2, mar = c(0, 0, 2, 0),
    text.col = c('black', venn_col[3:2], 'black', venn_col[c(1, 4)], 'black',
        'black')#, oma = rep(0, 4), pty = 'm'
)
dev.off()

## Venn by region width
venn <- vennCounts(countTable > 0)
print('annotation breakdown: all types')
cbind(venn, 'Percent' = round(venn[, 4]/ sum(venn[, 4]) * 100, 2))
venn[, 4] <- round(venn[, 4]/ sum(venn[, 4]) * 100, 2)

pdf(file.path(plotdir, paste0(opt$histone, '_ensemblVenn_percent.pdf')),
    width = 10, height = 10)
vennDiagram_custom(venn, 
    main = 'dbPeaks overlap with Ensembl v75 features (in %)', cex.main = 2,
    circle.col = venn_col[1:3], lwd = 1.5, cex = 2, mar = c(0, 0, 2, 0),
    text.col = c('black', venn_col[3:2], 'black', venn_col[c(1, 4)], 'black',
        'black')#, oma = rep(0, 4), pty = 'm'
)
dev.off()


getCount <- function(ann) {
    as.integer(countTable[, ann] > 0)
}
venn_idx <- lapply(c('exon', 'intergenic', 'intron'), getCount)
names(venn_idx) <- c('exon', 'intergenic', 'intron')
venn_bp <- venn
for(i in which(venn_bp[, 4] > 0)) {
    venn_bp[i, 4] <- round(sum(width(regions)[ venn_idx$exon == venn_bp[i, 1] & venn_idx$intergenic == venn_bp[i, 2] & venn_idx$intron == venn_bp[i, 3] ]) / 1e6, 2)
}

pdf(file.path(plotdir, paste0(opt$histone, '_ensemblVenn_bp.pdf')), width = 10,
    height = 10)
vennDiagram_custom(venn_bp, 
    main = 'dbPeaks overlap with Ensembl v75 features (in Mb)', cex.main = 2,
    circle.col = venn_col[1:3], lwd = 1.5, cex = 2, mar = c(0, 0, 2, 0),
    text.col = c('black', venn_col[3:2], 'black', venn_col[c(1, 4)], 'black',
        'black')#, oma = rep(0, 4), pty = 'm'
)
dev.off()
print('annotation breakdown: by region width')
cbind(venn_bp, 'Percent' = round(venn_bp[, 4]/ sum(venn_bp[, 4]) * 100, 2))
venn_bp[, 4] <- round(venn_bp[, 4]/ sum(venn_bp[, 4]) * 100, 2)

print('Percent fold change: number / bp')
cbind(venn_bp[, 1:3], 'Fold Change #/bp' = round(venn[, 4] / venn_bp[, 4], 2))
print('Percent fold change: bp / number')
cbind(venn_bp[, 1:3], 'Fold Change bp/#' = round(venn_bp[, 4] / venn[, 4], 2))



pdf(file.path(plotdir, paste0(opt$histone, '_ensemblVenn_bp_percent.pdf')),
    width = 10, height = 10)
vennDiagram_custom(venn_bp, 
    main = 'dbPeaks overlap with Ensembl v75 features (bp %)', cex.main = 2,
    circle.col = venn_col[1:3], lwd = 1.5, cex = 2, mar = c(0, 0, 2, 0),
    text.col = c('black', venn_col[3:2], 'black', venn_col[c(1, 4)], 'black',
        'black')#, oma = rep(0, 4), pty = 'm'
)
dev.off()

## Load coverage matrix
message(paste(Sys.time(), 'loading coverageMatrix.Rdata'))
load(file.path(maindir, 'coverageMatrix.Rdata'))

## Subset to regions of interest
coverageMatrix <- coverageMatrix[keepIndex, ]

message(paste(Sys.time(), 'log2 transforming counts'))
y <- log2(coverageMatrix + 1)

## PCA analysis
message(paste(Sys.time(), 'performing PCA analysis'))
pcList <- lapply(annoClassList, function(ii) {
	cat(".")
	pc = prcomp(t(y[ii,]))
	pc$rot = NULL # drop rotations
	return(pc) 
})

getPcaVars <- function(pca)  signif(((pca$sdev)^2)/(sum((pca$sdev)^2)),3)*100
pcVarMat <- sapply(pcList, getPcaVars)
rownames(pcVarMat) <- paste0("PC", seq_len(nrow(pcVarMat)))
pc1Mat <- sapply(pcList, function(x) x$x[,1])
pc2Mat <- sapply(pcList, function(x) x$x[,2])

## Make PCA plots
name <- c('Exonic', 'Intronic', 'Intergenic', 'Exon+Intron', 'All')
group <- factor(paste0(pd$BrainRegion, '_', c('52-',
    '52+')[as.numeric(pd$AgeDeath < 52) + 1]))
cellgroup <- factor(pd$CellType, levels = c('NeuN-', 'NeuN+'))
groupSimple <- groupInfo
levels(groupSimple) <- gsub('\\[23,52\\)', '52-', levels(groupSimple))
levels(groupSimple) <- gsub('\\[52,65\\]', '52+', levels(groupSimple))
levels(groupSimple) <- gsub('_', ':', levels(groupSimple))

pdf(file.path(plotdir, paste0(opt$histone, '_dbPeaks_PCA_byAnno.pdf')))
palette(brewer.pal(4, 'Paired'))
par(mar=c(5,6,2,2))
for(i in 1:ncol(pc1Mat)) {
	plot(x=pc1Mat[,i], y=pc2Mat[,i],
        bg = as.numeric(group),
		pch = c(21,22)[as.numeric(cellgroup)],
		xlab = paste0("PC1: ",pcVarMat[1,i],"% of Var Expl"),
		ylab = paste0("PC2: ",pcVarMat[2,i],"% of Var Expl"),
		cex.axis=2,cex.lab=2, cex.main=1.8,
		main = paste0("PCA of dbPeaks (", name[i],")"))
	legend("bottomright", c('Neun-', 'Neun+'), 
		pch=c(19,15), cex=1.5, ncol = 2, bty = 'n')
    legend("topright", levels(group),
        col = seq_len(length(levels(group))), lwd=5, cex=1.5, ncol = 2,
        bty = 'n')
}
dev.off()

pdf(file.path(plotdir, paste0(opt$histone, '_dbPeaks_PCsbyGroup_byAnno.pdf')),
    width = 11)
palette(brewer.pal(4, 'Paired'))
par(mar=c(14,6,2,2))
set.seed(20160514)
for(i in seq_len(ncol(pc1Mat))) {
	## PC1
	boxplot(pc1Mat[, i] ~ groupSimple, las=3,
		ylab = paste0("PC1: ", pcVarMat[1, i], "% of Var Expl"),
		cex.axis=1.7, cex.lab=2, cex.main=1.8, xlab="", outline=FALSE,
		main = paste0("PCA of dbPeaks (", name[i],")"))
	points(pc1Mat[, i] ~ jitter(as.numeric(groupSimple), amount=0.2),
		bg = as.numeric(group), cex=1.3,
		pch = c(21,22)[as.numeric(cellgroup)])
	# PC2 
	boxplot(pc2Mat[, i] ~ groupSimple, las=3,
		ylab = paste0("PC2: ", pcVarMat[2, i], "% of Var Expl"),
		cex.axis=1.7, cex.lab=2, xlab="", outline=FALSE)
	points(pc2Mat[, i] ~ jitter(as.numeric(groupSimple), amount=0.2),
		bg = as.numeric(group), cex=1.3,
		pch = c(21,22)[as.numeric(cellgroup)])
}
dev.off()


## Joint modeling
message(paste(Sys.time(), 'performing joint modeling'))
system.time( sumSqList <- parallel::mclapply(seq_len(nrow(y)), function(i) {
	if(i %% 10000 == 0) cat(".")
        t(anova(lm(y[i,] ~ BrainRegion + CellType + AgeDeath + Hemisphere + PMI + pH + Sex + Height + Weight + ChromatinAmount + AntibodyAmount + totalMapped + Individual_ID + FlowcellBatch + LibraryBatch, data=pd))[2])
}, mc.cores=8) )

ssOut <- do.call("rbind", sumSqList)
rownames(ssOut) <- NULL
bg <- matrix(rep(rowSums(ssOut), ncol(ssOut)), 
	ncol = ncol(ssOut), nrow = nrow(ssOut))
ssMat <- ssOut / bg
lab <- c('Brain region', 'Cell type', 'Age at death', 'Hemisphere', 'PMI', 'pH', 'Sex', 'Height', 'Weight', 'Chromatin amount', 'Antibody amount', 'Mapped reads', 'Individual', 'Flowcell batch', 'Library batch', 'Residual variation')

message(paste(Sys.time(), 'saving joint modeling results'))
save(ssMat, lab, file = file.path(maindir, paste0('ssMat_', opt$histone,
    '.Rdata')), compress=TRUE)

## Boxplot by variables
pdf(file.path(plotdir, paste0(opt$histone, '_boxplots_overall.pdf')),
    width = 12, height = 5)
par(mar=c(9,5,2,2))
palette(brewer.pal(7, "Dark2"))
boxplot(100*ssMat, xaxt="n", ylim = c(0, 100), 
	cex.axis=1.3,cex.lab=1.1, range=2,
	ylab="Percentage variance explained", cex=0.5)
text(seq_len(ncol(ssMat)) + 0.2, y = -8, lab, xpd=TRUE, srt=45, pos=2)
text(x = 8.5, y= 90, "All Regions", cex=1.7)
for(i in seq(along = annoClassList[1:4])) {
	ii = annoClassList[[i]]
	boxplot(100 * ssMat[ii,],xaxt="n", ylim = c(0, 100), 
		cex.axis=1.3, cex.lab=1.1, range=2, col = i,
		ylab="Percentage variance explained", cex=0.5)
	text(seq_len(ncol(ssMat)) + 0.1, y = -8, lab, xpd=TRUE, srt=45, pos=2)
	text(x = 8.5, y= 90, name[i], cex=1.7)
}
dev.off()

## Wihtout CellType
message(paste(Sys.time(), 'performing joint modeling without CellType'))
system.time( sumSqList2 <- parallel::mclapply(seq_len(nrow(y)), function(i) {
	if(i %% 10000 == 0) cat(".")
        t(anova(lm(y[i,] ~ BrainRegion + AgeDeath + Hemisphere + PMI + pH + Sex + Height + Weight + ChromatinAmount + AntibodyAmount + totalMapped + Individual_ID + FlowcellBatch + LibraryBatch, data=pd))[2])
}, mc.cores=8) )

ssOut2 <- do.call("rbind", sumSqList2)
rownames(ssOut2) <- NULL
bg2 <- matrix(rep(rowSums(ssOut2), ncol(ssOut2)), 
	ncol = ncol(ssOut2), nrow = nrow(ssOut2))
ssMat2 <- ssOut2 / bg2
lab2 <- c('Brain region', 'Age at death', 'Hemisphere', 'PMI', 'pH', 'Sex', 'Height', 'Weight', 'Chromatin amount', 'Antibody amount', 'Mapped reads', 'Individual', 'Flowcell batch', 'Library batch', 'Residual variation')

message(paste(Sys.time(), 'saving joint modeling results'))
save(ssMat2, lab2, file = file.path(maindir, paste0('ssMat2_', opt$histone,
    '_noCellType.Rdata')), compress=TRUE)

## Boxplot by variables
pdf(file.path(plotdir, paste0(opt$histone, '_boxplots_overall_noCellType.pdf')),
    width = 12, height = 5)
par(mar=c(9,5,2,2))
palette(brewer.pal(7, "Dark2"))
boxplot(100*ssMat2, xaxt="n", ylim = c(0, 100), 
	cex.axis=1.3,cex.lab=1.1, range=2,
	ylab="Percentage variance explained", cex=0.5)
text(seq_len(ncol(ssMat2)) + 0.2, y = -8, lab2, xpd=TRUE, srt=45, pos=2)
text(x = 8.5, y= 90, "All Regions", cex=1.7)
for(i in seq(along = annoClassList[1:4])) {
	ii = annoClassList[[i]]
	boxplot(100 * ssMat2[ii,],xaxt="n", ylim = c(0, 100), 
		cex.axis=1.3, cex.lab=1.1, range=2, col = i,
		ylab="Percentage variance explained", cex=0.5)
	text(seq_len(ncol(ssMat2)) + 0.1, y = -8, lab2, xpd=TRUE, srt=45, pos=2)
	text(x = 8.5, y= 90, name[i], cex=1.7)
}
dev.off()

## Load mean coverage
message(paste(Sys.time(), 'loading meanCoverage.Rdata'))
load(file.path(maindir, 'meanCoverage.Rdata'))

message(paste(Sys.time(), 'subsetting meanCoverage'))
map <- seq_len(length(fullRegions))[keepIndex]
meanCoverage <- meanCoverage[as.character(map)]
names(meanCoverage) <- seq_len(length(meanCoverage))


## Select regions to highlight per covariate
highlight <- apply(ssMat, 2, function(x) {
    head(order(x, decreasing = TRUE), 10)
})


## From https://github.com/leekgroup/derSupplement/commit/d67dc1d2aee34eaae6e4a63082d03ea4397d46f1
## Load info
sql_file <- "/home/epi/ajaffe/Lieber/Projects/RNAseq/Ribozero_Compare/TxDb.Hsapiens.BioMart.ensembl.GRCh37.p12/inst/extdata/TxDb.Hsapiens.BioMart.ensembl.GRCh37.p12.sqlite"
TranscriptDb <- loadDb(sql_file)

## Fix seqlevels
seqlevels(TranscriptDb, force=TRUE) <- c(1:22,"X","Y","MT")
seqlevels(TranscriptDb) <- paste0("chr", c(1:22,"X","Y","M"))
ensGene <- genes(TranscriptDb)

## Get data from BioMart
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", # VERSION 75, hg19
 	dataset="hsapiens_gene_ensembl",
	host="feb2014.archive.ensembl.org")
sym <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
 	values = names(ensGene), mart = ensembl)
ensGene$Symbol <- sym$hgnc_symbol[match(names(ensGene), sym$ensembl_gene_id)]
ensGene <- ensGene[!grepl("^MIR[0-9]", ensGene$Symbol)] # drop mirs


## Get nearest annotation
genes <- annotateTranscripts(txdb = TranscriptDb)
matchGenes(x = zoom, subject = genes)

shorten <- function(x) {
    x <- gsub('ACC', 'A', x)
    x <- gsub('DLPFC', 'D', x)
    x <- gsub('52', '', x)
    x <- gsub('NeuN', 'N', x)
    return(x)
}

groupSimple_mean <- factor(shorten(levels(groupSimple)), levels = shorten(levels(groupSimple)))

message(paste(Sys.time(), 'creating highlight region plots'))
for(h in seq_len(ncol(highlight))) {
    message(paste(Sys.time(), 'highlighting', colnames(ssMat)[h]))
    ## Subset data
    regs <- regions[highlight[, h]]
    meanC <- meanCoverage[as.character(highlight[, h])]
    names(meanC) <- seq_len(10)
    annRegs <- list(annotationList = ensemblAnno$annotationList[highlight[, h]])
    names(annRegs$annotationList) <- seq_len(10)
    for(j in seq_len(10)) {
        ## Remove symbols
        annRegs[[1]][[j]]$symbol <- CharacterList('')
    }    
    nearestAnn <- matchGenes(x = regs, subject = genes)
    nearestAnn$name <- NA
    ov <- findOverlaps(ensGene, regs, ignore.strand = TRUE)
    nearestAnn$name[subjectHits(ov)] <- ensGene$Symbol[queryHits(ov)]
    
    pdf(file.path(plotdir, paste0(opt$histone, '_highlight_',
        colnames(ssMat)[h], '.pdf')), height = 5, width = 8)
    plotRegionCoverage(regions = regs,
        regionCoverage = meanC,
        groupInfo = groupSimple_mean,
        nearestAnnotation = nearestAnn,
        annotatedRegions = annRegs,
        ask = FALSE, verbose = FALSE,
        txdb = TranscriptDb, colors = brewer.pal(12, 'Paired')[5:12]
    )
    dev.off()
}



## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
