## Analyze results

## Load libraries
library('getopt')

## Available at http://www.bioconductor.org/packages/release/bioc/html/derfinder.html
library('derfinder')
library('GenomicRanges')
library('limma')
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
plotdir <- file.path(maindir, 'plots')
dir.create(plotdir, showWarnings = FALSE)
    
## Results
load(file.path(maindir, 'groupInfo.Rdata'))
load(file.path(maindir, 'colsubset.Rdata'))
## Load regions
message(paste(Sys.time(), 'loading fullRegions.Rdata'))
load(file.path(maindir, 'fullRegions.Rdata'))

keepIndex <- width(fullRegions) >=6
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

pdf(file.path(plotdir, 'ensemblVenn.pdf'), width = 10, height = 10)
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

pdf(file.path(plotdir, 'ensemblVenn_percent.pdf'), width = 10, height = 10)
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
for(i in which(venn[, 4] > 0)) {
    venn[i, 4] <- round(sum(width(regions)[ venn_idx$exon == venn[i, 1] & venn_idx$intergenic == venn[i, 2] & venn_idx$intron == venn[i, 3] ]) / 1e6, 2)
}

pdf(file.path(plotdir, 'ensemblVenn_bp.pdf'), width = 10, height = 10)
vennDiagram_custom(venn, 
    main = 'dbPeaks overlap with Ensembl v75 features (in Mb)', cex.main = 2,
    circle.col = venn_col[1:3], lwd = 1.5, cex = 2, mar = c(0, 0, 2, 0),
    text.col = c('black', venn_col[3:2], 'black', venn_col[c(1, 4)], 'black',
        'black')#, oma = rep(0, 4), pty = 'm'
)
dev.off()
print('annotation breakdown: by region width')
cbind(venn, 'Percent' = round(venn[, 4]/ sum(venn[, 4]) * 100, 2))

venn[, 4] <- round(venn[, 4]/ sum(venn[, 4]) * 100, 2)

pdf(file.path(plotdir, 'ensemblVenn_bp_percent.pdf'), width = 10, height = 10)
vennDiagram_custom(venn, 
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
pdf(file.path(plotdir, 'dbPeaks_PCA_byAnno.pdf'))
palette(brewer.pal(4, "Dark2"))
name <- c("Exonic", "Intronic", "Intergenic","Exon+Intron", "All")
group <- factor(paste0(pd$BrainRegion, '_', c('52-',
    '52+')[as.numeric(pd$AgeDeath < 52) + 1]))
par(mar=c(5,6,2,2))
for(i in 1:ncol(pc1Mat)) {
	plot(x=pc1Mat[,i], y=pc2Mat[,i],
        bg = as.numeric(group),
		pch = c(21,22)[as.numeric(factor(pd$CellType, levels = c('NeuN-', 'NeuN+')))],
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




## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
