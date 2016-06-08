## Analyze results

options(width = 180)

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
library('Hmisc')
library('GOstats')
library('org.Hs.eg.db')
library('devtools')
library('ggplot2')

## Specify parameters
spec <- matrix(c(
    'histone', 'i', 1, 'character', 'For epimap, the histone mark to use. Either H3K27ac or H3K4me3.',
    'cutoff', 't', 1, 'numeric', 'Cutoff to use',
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
if(FALSE) {
    opt <- list(histone = 'H3K4me3', cutoff = 10)
    opt <- list(histone = 'H3K27ac', cutoff = 10)
}
    
## Number of cores used
cores <- 10

## Dirs
rootdir <- '/dcl01/lieber/ajaffe/derRuns/derChIP/epimap/'
maindir <- file.path(rootdir, 'derAnalysis', paste0('run1-v1.5.38-',
    opt$histone))
regmatdir <- file.path(rootdir, 'regionMatrix')
resdir <- file.path(regmatdir, paste0(opt$histone, '-cut', opt$cutoff))
dir.create(resdir, showWarnings = FALSE)
plotdir <- file.path(regmatdir, paste0('plots-cut', opt$cutoff))
dir.create(plotdir, showWarnings = FALSE)

## Age cutoff
ageCut <- ifelse(opt$histone == 'H3K4me3', 52, 53)

## Results
load(file.path(maindir, 'groupInfo.Rdata'))
load(file.path(maindir, 'colsubset.Rdata'))
load(file.path(maindir, 'models.Rdata'))

## Load regions
message(paste(Sys.time(), 'loading fullRegions.Rdata'))
load(file.path(regmatdir, paste0('fullRegions-', opt$histone, '-cut',
    opt$cutoff, '.Rdata')))

## Load coverage matrix
message(paste(Sys.time(), 'loading coverageMatrix.Rdata'))
load(file.path(regmatdir, paste0('coverageMatrix-', opt$histone, '-cut',
    opt$cutoff, '.Rdata')))

message(paste(Sys.time(), 'log2 transforming counts'))
y <- log2(coverageMatrix + 1)

# get the f statistic from 2 lmFit objects
getF <- function(fit, fit0, theData) {
	
	rss1 = rowSums((fitted(fit)-theData)^2)
	df1 = ncol(fit$coef)
	rss0 = rowSums((fitted(fit0)-theData)^2)
	df0 = ncol(fit0$coef)

	fstat = ((rss0-rss1)/(df1-df0))/(rss1/(ncol(theData)-df1))
	f_pval = pf(fstat, df1-1, ncol(theData)-df1,lower.tail=FALSE)
	fout = cbind(fstat,df1-1,ncol(theData)-df1,f_pval)
	colnames(fout)[2:3] = c("df1","df0")
	fout = data.frame(fout)
	return(fout)
}

## Determine which ERs are DERs
message(paste(Sys.time(), 'performing DE analysis'))
fit <- lmFit(y, models$mod)
eb <- ebayes(fit)
fit0 <- lmFit(y, models$mod0)
ff <- getF(fit, fit0, y)

## Adjust p-values by Bonf
message(paste(Sys.time(), 'adjusting p-values'))
ff$p_bonf <- p.adjust(ff$f_pval, "bonferroni")
ff$significantFWER <- ifelse(ff$p_bonf < 0.05, 'TRUE', 'FALSE')
mcols(fullRegions) <- cbind(mcols(fullRegions), ff)

message(paste(Sys.time(), 'saving new fullRegions object'))
save(fullRegions, file = file.path(resdir, 'fullRegions.Rdata'))

message(paste(Sys.time(), 'subsetting the regions and coverage matrix'))
#keepIndex <- width(fullRegions) >= 6
keepIndex <- fullRegions$significantFWER == 'TRUE'
regions <- fullRegions[keepIndex]
print('Number of ERs, DERs and percent of ERs that are DERs')
length(fullRegions)
length(regions)
round(length(regions) / length(fullRegions) * 100, 2)

## Subset to regions of interest
coverageMatrix <- coverageMatrix[keepIndex, ]
y <- y[keepIndex, ]

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
save(ensemblAnno, file = file.path(resdir, 'ensemblAnno.Rdata'))
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
group <- factor(paste0(pd$BrainRegion, ':', paste0(ageCut, c('-', '+'))[as.numeric(pd$AgeDeath < ageCut) + 1]))
cellgroup <- factor(pd$CellType, levels = c('NeuN-', 'NeuN+'))
groupSimple <- groupInfo
levels(groupSimple) <- gsub(paste0('\\[23,', ageCut, '\\)'), paste0(ageCut, '-'), levels(groupSimple))
levels(groupSimple) <- gsub(paste0('\\[', ageCut, ',65\\]'), paste0(ageCut, '+'), levels(groupSimple))
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
		main = paste0("PCA of dbPeaks (", name[i],")"),
        ylim = range(pc2Mat[, i]) * 1.2)
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

## Add BMI
pd$BMI <- pd$Weight / (pd$Height^2) * 703
print('Number of samples, height, weight and BMI by sex')
table(pd$Sex)
tapply(pd$Height, pd$Sex, summary)
tapply(pd$Weight, pd$Sex, summary)
tapply(pd$BMI, pd$Sex, summary)

## Joint modeling
message(paste(Sys.time(), 'performing joint modeling'))
system.time( sumSqList <- parallel::mclapply(seq_len(nrow(y)), function(i) {
	if(i %% 10000 == 0) cat(".")
        t(anova(lm(y[i,] ~ BrainRegion + CellType + AgeDeath + Hemisphere + PMI + pH + Sex + Height + BMI + ChromatinAmount + totalMapped + Individual_ID + FlowcellBatch + LibraryBatch, data=pd))[2])
}, mc.cores = cores) )

ssOut <- do.call("rbind", sumSqList)
rownames(ssOut) <- NULL
bg <- matrix(rep(rowSums(ssOut), ncol(ssOut)), 
	ncol = ncol(ssOut), nrow = nrow(ssOut))
ssMat <- ssOut / bg
lab <- c('Brain region', 'Cell type', 'Age at death', 'Hemisphere', 'PMI', 'pH', 'Sex', 'Height', 'BMI', 'Chromatin amount', 'Mapped reads', 'Individual', 'Flowcell batch', 'Library batch', 'Residual variation')
names(lab) <- colnames(ssMat)

message(paste(Sys.time(), 'saving joint modeling results'))
save(ssMat, lab, file = file.path(resdir, paste0('ssMat_', opt$histone,
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

## Without CellType
message(paste(Sys.time(), 'performing joint modeling without CellType'))
system.time( sumSqList2 <- parallel::mclapply(seq_len(nrow(y)), function(i) {
	if(i %% 10000 == 0) cat(".")
        t(anova(lm(y[i,] ~ BrainRegion + AgeDeath + Hemisphere + PMI + pH + Sex + Height + BMI + ChromatinAmount + totalMapped + Individual_ID + FlowcellBatch + LibraryBatch, data=pd))[2])
}, mc.cores = cores) )

ssOut2 <- do.call("rbind", sumSqList2)
rownames(ssOut2) <- NULL
bg2 <- matrix(rep(rowSums(ssOut2), ncol(ssOut2)), 
	ncol = ncol(ssOut2), nrow = nrow(ssOut2))
ssMat2 <- ssOut2 / bg2
lab2 <- c('Brain region', 'Age at death', 'Hemisphere', 'PMI', 'pH', 'Sex', 'Height', 'BMI', 'Chromatin amount', 'Mapped reads', 'Individual', 'Flowcell batch', 'Library batch', 'Residual variation')

message(paste(Sys.time(), 'saving joint modeling results (no cellType)'))
save(ssMat2, lab2, file = file.path(resdir, paste0('ssMat2_', opt$histone,
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


## p-values: main 3 covariates
message(paste(Sys.time(), 'ANOVA main 3 covariates'))
i_groups <- split(seq_len(nrow(y)), Hmisc::cut2(seq_len(nrow(y)), m = 5e3))
system.time( pList <- parallel::mclapply(i_groups, function(i_group) {
    res <- data.frame(matrix(NA, ncol = 3, nrow = length(i_group)))
    colnames(res) <- c('BrainRegion', 'CellType', 'AgeDeath')
    
    for(j in seq_len(length(i_group))) {
        i <- i_group[j]
    	if(i %% 10000 == 0) cat(".")
           res[j, ] <- c(
                'BrainRegion' = anova(lm(y[i, ] ~ BrainRegion, data = pd))[, 'Pr(>F)'][1],
                'CellType' = anova(lm(y[i, ] ~ CellType, data = pd))[, 'Pr(>F)'][1],
                'AgeDeath' = anova(lm(y[i, ] ~ AgeDeath, data = pd))[, 'Pr(>F)'][1]
            )
    }
	return(res)
}, mc.cores = cores) )
system.time( pTable_main <- do.call(rbind, pList) )
rownames(pTable_main) <- NULL
message(paste(Sys.time(), 'saving pTable_main.Rdata'))
save(pTable_main, file = file.path(resdir, 'pTable_main.Rdata'))

## Summary info for p-value table
print('Summary pTable_main')
summary(pTable_main)
print('Summary pTable_main (Bonf adjusted)')
pTable_main_adj <- t(apply(pTable_main, 1, p.adjust, method = 'bonferroni'))
summary(pTable_main_adj)
print('Summary of -log10 pvalues (Bonferroni adjusted)')
summary(-log(pTable_main_adj))
pTable_main_sig <- pTable_main_adj < 0.05
print('Significant dbPeaks by modeled covariate')
colSums(pTable_main_sig)
round(colSums(pTable_main_sig) / nrow(pTable_main_sig) * 100, 2)

venn_main <- vennCounts(pTable_main_sig)
colnames(venn_main) <- c('Brain region', 'Cell type', 'Age at death', 'Counts')
print('dbPeaks by main covariate')
venn_main
venn_col <- brewer.pal(4, "Set1")[2:4]
pdf(file.path(plotdir, paste0(opt$histone, '_venn_mainCovariates.pdf')),
    width = 10, height = 10)
vennDiagram_custom(venn_main, 
    main = paste(opt$histone, 'dbPeaks by main covariates'), cex.main = 2,
    circle.col = venn_col[1:3], lwd = 1.5, cex = 2, mar = c(0, 0, 2, 0),
    text.col = c('black', venn_col[3:2], 'black', venn_col[1], 'black', 'black',
        'black')#, oma = rep(0, 4), pty = 'm'
)
dev.off()

## Subset PCAs
covClassList <- lapply(c('BrainRegion', 'CellType', 'AgeDeath'), function(covariate) {
    not <- !pTable_main_sig[, -which(colnames(pTable_main_sig) == covariate)]
    which(pTable_main_sig[, covariate] & apply(not, 1, all))
})
names(covClassList) <- c('BrainRegion', 'CellType', 'AgeDeath')
print('Percent of regions exclusive by main covariate')
round(sapply(covClassList, length) / nrow(y) * 100, 2)
covClassList <- covClassList[sapply(covClassList, length) > 1]


pcListCov <- lapply(covClassList, function(ii) {
	cat(".")
	pc = prcomp(t(y[ii,]))
	pc$rot = NULL # drop rotations
	return(pc) 
})

pcVarMatCov <- lapply(pcListCov, getPcaVars)
pcVarMatCov <- do.call(cbind, lapply(pcVarMatCov, head, n = min(sapply(pcVarMatCov, length))))
rownames(pcVarMatCov) <- paste0("PC", seq_len(nrow(pcVarMatCov)))
pc1MatCov <- sapply(pcListCov, function(x) x$x[,1])
pc2MatCov <- sapply(pcListCov, function(x) x$x[,2])

## Plots by covariate
pdf(file.path(plotdir, paste0(opt$histone, '_dbPeaks_PCA_byCovariate.pdf')))
palette(brewer.pal(4, 'Paired'))
par(mar=c(5,6,2,2))
for(i in seq_len(ncol(pc1MatCov))) {
	plot(x=pc1MatCov[,i], y=pc2MatCov[,i],
        bg = as.numeric(group),
		pch = c(21,22)[as.numeric(cellgroup)],
		xlab = paste0("PC1: ",pcVarMatCov[1,i],"% of Var Expl"),
		ylab = paste0("PC2: ",pcVarMatCov[2,i],"% of Var Expl"),
		cex.axis=2, cex.lab=2, cex.main=1.8,
		main = paste0("PCA of dbPeaks (", colnames(venn_main)[i],")"),
        ylim = range(pc2MatCov[,i]) * 1.4)
	legend("bottomright", c('Neun-', 'Neun+'), 
		pch=c(19,15), cex=1.5, ncol = 2, bty = 'n')
    legend("topright", levels(group),
        col = seq_len(length(levels(group))), lwd=5, cex=1.5, ncol = 2,
        bty = 'n')
}
dev.off()

## Plots by group and covariate
pdf(file.path(plotdir, paste0(opt$histone,
    '_dbPeaks_PCsbyGroup_byCovariate.pdf')), width = 11)
palette(brewer.pal(4, 'Paired'))
par(mar=c(14,6,2,2))
set.seed(20160516)
for(i in seq_len(ncol(pc1MatCov))) {
	## PC1
	boxplot(pc1MatCov[, i] ~ groupSimple, las=3,
		ylab = paste0("PC1: ", pcVarMatCov[1, i], "% of Var Expl"),
		cex.axis=1.7, cex.lab=2, cex.main=1.8, xlab="", outline=FALSE,
		main = colnames(venn_main)[i],
        ylim = range(pc1MatCov[, i]) * 1.1)
	points(pc1MatCov[, i] ~ jitter(as.numeric(groupSimple), amount=0.2),
		bg = as.numeric(group), cex=1.3,
		pch = c(21,22)[as.numeric(cellgroup)])
	# PC2 
	boxplot(pc2MatCov[, i] ~ groupSimple, las=3,
		ylab = paste0("PC2: ", pcVarMatCov[2, i], "% of Var Expl"),
		cex.axis=1.7, cex.lab=2, xlab="", outline=FALSE)
	points(pc2MatCov[, i] ~ jitter(as.numeric(groupSimple), amount=0.2),
		bg = as.numeric(group), cex=1.3,
		pch = c(21,22)[as.numeric(cellgroup)])
}
dev.off()


## p-values: remaining covariates
message(paste(Sys.time(), 'ANOVA remaining covariates'))
system.time( pList2 <- parallel::mclapply(i_groups, function(i_group) {
    res <- data.frame(matrix(NA, ncol = length(lab) - 4, nrow = length(i_group)))
    colnames(res) <- c('Hemisphere', 'PMI', 'pH', 'Sex', 'Height', 'BMI', 'ChromatinAmount', 'totalMapped', 'Individual_ID', 'FlowcellBatch', 'LibraryBatch')
    
    for(j in seq_len(length(i_group))) {
        i <- i_group[j]
    
    	if(i %% 10000 == 0) cat(".")
            fit1 <- lm(y[i, ] ~ BrainRegion + CellType + AgeDeath, data = pd)
            res[j, ] <- c(
                'Hemisphere' = anova(fit1, lm(y[i, ] ~ BrainRegion + CellType + AgeDeath + Hemisphere, data = pd))[, 'Pr(>F)'][2],
                'PMI' = anova(fit1, lm(y[i, ] ~ BrainRegion + CellType + AgeDeath + PMI, data = pd))[, 'Pr(>F)'][2],
                'pH' = anova(fit1, lm(y[i, ] ~ BrainRegion + CellType + AgeDeath + pH, data = pd))[, 'Pr(>F)'][2],
                'Sex' = anova(fit1, lm(y[i, ] ~ BrainRegion + CellType + AgeDeath + Sex, data = pd))[, 'Pr(>F)'][2],
                'Height' = anova(fit1, lm(y[i, ] ~ BrainRegion + CellType + AgeDeath + Height, data = pd))[, 'Pr(>F)'][2],
                'BMI' = anova(fit1, lm(y[i, ] ~ BrainRegion + CellType + AgeDeath + BMI, data = pd))[, 'Pr(>F)'][2],
                'ChromatinAmount' = anova(fit1, lm(y[i, ] ~ BrainRegion + CellType + AgeDeath + ChromatinAmount, data = pd))[, 'Pr(>F)'][2],
                'totalMapped' = anova(fit1, lm(y[i, ] ~ BrainRegion + CellType + AgeDeath + totalMapped, data = pd))[, 'Pr(>F)'][2],
                'Individual_ID' = anova(fit1, lm(y[i, ] ~ BrainRegion + CellType + AgeDeath + Individual_ID, data = pd))[, 'Pr(>F)'][2],
                'FlowcellBatch' = anova(fit1, lm(y[i, ] ~ BrainRegion + CellType + AgeDeath + FlowcellBatch, data = pd))[, 'Pr(>F)'][2],
                'LibraryBatch' = anova(fit1, lm(y[i, ] ~ BrainRegion + CellType + AgeDeath + LibraryBatch, data = pd))[, 'Pr(>F)'][2]
            )
    }
    return(res)
}, mc.cores = cores) )
system.time( pTable_rest <- do.call(rbind, pList2) )
rownames(pTable_rest) <- NULL
message(paste(Sys.time(), 'saving pTable_rest.Rdata'))
save(pTable_rest, file = file.path(resdir, 'pTable_rest.Rdata'))

## Summary info for p-value table
print('Summary pTable_rest')
summary(pTable_rest)
print('Summary pTable_rest (Bonf adjusted)')
pTable_rest_adj <- t(apply(pTable_rest, 1, p.adjust, method = 'bonferroni'))
summary(pTable_rest_adj)
print('Summary of -log10 pvalues (Bonferroni adjusted)')
summary(-log(pTable_rest_adj))
pTable_rest_sig <- pTable_rest_adj < 0.05
print('Significant dbPeaks by un-modeled covariate')
colSums(pTable_rest_sig)
round(colSums(pTable_rest_sig) / nrow(pTable_rest_sig) * 100, 2)



## Cluster un-modeled covariates
pdf(file.path(plotdir, paste0(opt$histone, '_clus_otherCovariates.pdf')))
plot(hclust(as.dist(1 - cor( -log(pTable_rest_adj) ))), frame.plot = TRUE, main = 'Un-modeled covariates clustered by -log10 p-value', ylab = '', axes = FALSE, sub = '', ann = TRUE, labels = lab[4:(length(lab) - 1)], xlab = '')
dev.off()


## Groups
venn_rest <- vennCounts(pTable_rest_sig)
counts <- venn_rest[, ncol(pTable_rest) + 1]
venn_not0 <- venn_rest[as.integer(names(sort(counts[counts > 0], decreasing = TRUE))), ]
head(venn_not0, 10)

## Define region sets
regSets <- apply(cbind(pTable_main_sig, pTable_rest_sig), 2, which)
regSets <- regSets[sapply(regSets, length) > 0]
message(paste(Sys.time(), 'saving regSets.Rdata'))
save(regSets, file = file.path(resdir, 'regSets.Rdata'))

print('Number of regions per covariate set')
sets_l <- sapply(regSets, length)
sets_l

## Repeated ones by set
repeated_sets <- sapply(regSets, function(s) {
    sapply(regSets, function(r) {
        sum(s %in% r)
    })
})
print('Regions that are repeated in other sets')
diag(repeated_sets) <- NA
repeated_sets
repeated_sets > 0
print('Regions that are repeated in other sets: percent by row')
round(repeated_sets / matrix(rep(sets_l, ncol(repeated_sets)), ncol = ncol(repeated_sets)) * 100, 2)
print('Regions that are repeated in other sets: percent by column')
round(repeated_sets / matrix(rep(sets_l, nrow(repeated_sets)), nrow = nrow(repeated_sets), byrow = TRUE) * 100, 2)



## Select regions to highlight per covariate
pTable <- -log10(cbind(pTable_main_adj, pTable_rest_adj))
highlight <- lapply(names(regSets), function(x) {
    ord <- order(pTable[, x][regSets[[x]]], decreasing = TRUE)
    ## Maximum 100
    res <- regSets[[x]][head(ord, 50)] # Limiting to 50 in case some regions are huge
    return(res)
})
names(highlight) <- names(regSets)
message(paste(Sys.time(), 'saving highlight.Rdata'))
save(highlight, file = file.path(resdir, 'highlight.Rdata'))

print('Number of regions highlighted by covariate')
high_l <- sapply(highlight, length)
high_l

## Repeated ones?
repeated <- sapply(highlight, function(s) {
    sapply(highlight, function(r) {
        sum(s %in% r)
    })
})
print('Highlighted regions that are repeated in other sets')
diag(repeated) <- NA
repeated
repeated > 0
print('Highlighted regions that are repeated in other sets: percent by row')
round(repeated / matrix(rep(high_l, ncol(repeated)), ncol = ncol(repeated)) * 100, 2)
print('Highlighted regions that are repeated in other sets: percent by column')
round(repeated / matrix(rep(high_l, nrow(repeated)), nrow = nrow(repeated), byrow = TRUE) * 100, 2)

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

shorten <- function(x) {
    x <- gsub('ACC', 'A', x)
    x <- gsub('DLPFC', 'D', x)
    x <- gsub(ageCut, '', x)
    x <- gsub('NeuN', 'N', x)
    return(x)
}

groupSimple_mean <- factor(shorten(levels(groupSimple)), levels = shorten(levels(groupSimple)))
tIndexes <- split(seq_len(length(groupInfo)), groupInfo)
files <- pd$bamFile
names(files) <- pd$Sample_ID

message(paste(Sys.time(), 'creating highlight region plots'))
for(h in seq_len(length(highlight))) {
    message(paste(Sys.time(), 'highlighting', names(highlight)[h]))
    ## Subset data
    regs <- regions[highlight[[h]]]
    ## Add 100 bp padding to each side
    regs <- resize(regs, width(regs) + 200, fix = 'center')
    
    message(paste(Sys.time(), 'calculating mean coverage'))
    regionCov <- getRegionCoverage(regions = regs, files = files,
        totalMapped = pd$totalMapped, verbose = FALSE)
    meanC <- lapply(regionCov, function(x) {
    	sapply(tIndexes, function(ii) rowMeans(x[, ii]))
    })
    
    h_len <- length(highlight[[h]])
    names(meanC) <- seq_len(h_len)
    annRegs <- list(annotationList = ensemblAnno$annotationList[highlight[[h]]])
    names(annRegs$annotationList) <- seq_len(h_len)
    for(j in seq_len(h_len)) {
        ## Remove symbols
        annRegs[[1]][[j]]$symbol <- CharacterList('')
    }    
    nearestAnn <- matchGenes(x = regs, subject = genes)
    nearestAnn$name <- NA
    ov <- findOverlaps(ensGene, regs, ignore.strand = TRUE)
    nearestAnn$name[subjectHits(ov)] <- ensGene$Symbol[queryHits(ov)]
    nearestAnn$name[nearestAnn$name == ''] <- NA
    
    pdf(file.path(plotdir, paste0(opt$histone, '_highlight_',
        names(highlight)[h], '.pdf')), height = 5, width = 8)
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


message(paste(Sys.time(), 'creating highlight region plots (scatterplot)'))
for(h in seq_len(length(highlight))) {
    message(paste(Sys.time(), 'highlighting', names(highlight)[h], '(scatterplots)'))
    height <- ifelse(names(highlight[h]) == 'FlowcellBatch', 10, 7)
    pdf(file.path(plotdir, paste0(opt$histone, '_highlight_',
        names(highlight)[h], '_scatter.pdf')), height = height)
    for(hh in highlight[[h]]) {
        df <- data.frame(x = pd[, names(highlight)[h]], y = y[hh, ], group = group, cell = cellgroup)
        ypos <- max(df$y)
        xpos <- if(is.character(df$x) | is.factor(df$x)) length(unique(df$x)) - 0.5 else max(mean(range(df$x)), 0.9 * max(df$x))
        if(names(highlight)[h] == 'LibraryBatch') {
            xpos <- xpos * 0.9
            df$x <- factor(df$x, levels = sort(as.integer(levels(df$x))))
        }
        if(names(highlight)[h] == 'FlowcellBatch') xpos <- xpos - 1
        if(names(highlight)[h] == 'Individual_ID') xpos <- xpos - 2
        label <- paste('-log10 p-value', round(pTable[, names(highlight)[h]][hh], 2))
        g <- ggplot(df, aes(x = x, y = y, colour = group, shape = cell)) 
        g <- if(is.character(df$x) | is.factor(df$x)) g + geom_jitter(size = 3, width = 0.2, height = NULL) else g + geom_point(size = 3)
        g <- g +
            scale_color_manual(values = brewer.pal(4, 'Paired')) +
            scale_shape_manual(values = c(16, 15)) +
            theme_bw(base_size = 16) +
            xlab(lab[names(highlight)[h]]) +
            ylab('Counts (log2)') +
            annotate('text', label = label, y = ypos, x = xpos)
        g <- if((is.character(df$x) | is.factor(df$x)) & length(unique(df$x)) > 4) g + theme(legend.position = 'none', axis.text.x = element_text(angle = 90, hjust = 1)) else g + theme(legend.position = 'none')
        print(g)
    }
    dev.off()
}


## GO analysis
# Modified from /home/epi/ajaffe/Lieber/lieber_functions_aj.R
dogo <- function(names, Universe, goP = 0.01, cond = FALSE, ontology = 'BP'){
    gomap_names <- org.Hs.egREFSEQ2EG
    x <- unlist(mget(as.character(names), gomap_names, ifnotfound = NA))
    x <- x[!is.na(x)]    
    Universe <- unique(c(Universe, unique(x)))

    params <- new("GOHyperGParams", geneIds = unique(x),
                  universeGeneIds = Universe,
                  annotation = 'org.Hs.eg.db',
                  ontology = ontology, pvalueCutoff = goP, conditional = cond,
                  testDirection="over")
    ht <- hyperGTest(params)
    tab <- summary(ht)
    tmp1 <- geneIdsByCategory(ht)
    tmp1 <- tmp1[tab[, 1]]
    tab$IDs <- sapply(tmp1, function(y) paste(names(x)[x %in% y], collapse=";"))
    return(tab)
}

## Define universe: genes with a dbPeak within 5kb
bg_genes <- resize(ensGene, width(ensGene) + 1e4, fix = 'center')
bg_genes <- bg_genes[countOverlaps(bg_genes, regions, ignore.strand = TRUE) > 0]
print('Percent of genes included in the background')
round(length(bg_genes) / length(ensGene) * 100, 2)

## Define GO universe
gomap <- org.Hs.egENSEMBL2EG
bg_universe <- unlist(mget(as.character(bg_genes$gene_id), gomap, ifnotfound = NA))
print('Percent of background genes missing')
round(sum(is.na(bg_universe)) / length(bg_universe) * 100, 2)
bg_universe <- bg_universe[!is.na(bg_universe)]


genes.knownGene <- annotateTranscripts(txdb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)

message(paste(Sys.time(), 'performing GO analysis'))
goByCovariate <- lapply(regSets, function(ii) {
    regs <- regions[ii]
    annotation.knownGene <- matchGenes(x = regs, subject = genes.knownGene)
    regs_names <- unlist(strsplit(annotation.knownGene$annotation, ' '))
    ## Clean up
    regs_names <- regs_names[!is.na(regs_names)]
    go <- tryCatch(dogo(regs_names, bg_universe), error = function(e) return(NULL))
})
names(goByCovariate) <- names(regSets)

goByCovariate <- goByCovariate[!sapply(goByCovariate, is.null)]

message(paste(Sys.time(), 'saving goByCovariate.Rdata'))
save(goByCovariate, file = file.path(resdir, 'goByCovariate.Rdata'))

## Show GO results
print('Top GO results')
for(g in names(regSets)) {
    print(g)
    print(head(goByCovariate[[g]][, -8], 20))
}


## Reproducibility info
proc.time()
message(Sys.time())
devtools::session_info()
