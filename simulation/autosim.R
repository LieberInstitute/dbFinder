## The original version of this file is available from 
## http://bioinf.wehi.edu.au/csaw/
## and the simulation results are described in
## csaw: a Bioconductor package for differential binding analysis of ChIP-seq 
## data using sliding windows by Aaron Lun and Gordon K. Smyth
## Nucleic Acids Research 44, e45 (2015).
##
## Summary of modifications:
## 

# This simulation generates a mixture of complex peak structures.
# This is done by treating each peak as a mixture of three subpeaks.
# For DB peaks, one or more of those subpeaks are randomly chosen and eliminated.

library('csaw')
library('edgeR')
library('DiffBind')
library('derfinder')
library('devtools')
library('Rsamtools')
source('xscss_modified.R')

## For passing arguments
library('getopt')

## Specify parameters
spec <- matrix(c(
	'type', 't', 1, 'character', 'Either tf or hist',
    'clean', 'c', 2, 'logical', 'Whether to clean up',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## Default values
if (is.null(opt$clean)) opt$clean <- FALSE

if(FALSE) {
    ## For testing
    opt <- list(type = 'hist', clean = FALSE)
    it <- 1
    bam.files <- dir(pattern = 'bam$')
}

stopifnot(opt$type %in% c('tf', 'hist'))
is.tf <- opt$type == 'tf'

fraglen <- 100
iters <- 1
npeaks <- 20000
nde <- 500

prior.df <- 20
dispersion <- 1 / prior.df # For the sake of possible flexibility later; right now, it just cancels out.
grouping <- rep(c('A', 'B'), each = 5)
true.widths <- c(300, 500, 700)
width.n <- length(true.widths)
base.mu <- 30

if (is.tf) { 
	radius <- fraglen
	up.mu <- 45
	down.mu <- 15
} else {
	
    radiuses <- true.widths / 2L
}

################################################################################

design <- model.matrix(~factor(grouping))
fdr.thresholds <- c(0.01, 0.05, 0.1, 0.15, 0.2)

result.file <- paste0(all.fix, '_result.tsv')
unlink(result.file)
dump2file <- function(id, cutoff, result) {
	write.table(file = result.file, 
        data.frame(id, cutoff, 1 - result$overlap/result$found, result$recall), 
		sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE,
        append = TRUE)
}
set.seed(20160927)

################################################################################

for (it in seq_len(iters)) {
    ofile <- paste0(all.fix, '_it', it, '_out.tsv')
    
	### Generating simulated data for histone mark data.
	up.pk <- seq_len(nde)
	down.pk <- seq_len(nde) + nde
	disp <- prior.df * dispersion / rchisq(npeaks * width.n, df=prior.df)
    disp <- split(disp, rep(seq_len(width.n), each = npeaks))
    names(disp) <- true.widths


	type.A <- type.B <- lapply(seq_len(width.n), function(x) {
        lapply(seq_len(3), function(y) { !logical(npeaks) })})
    pos <- lapply(seq_len(width.n), function(x) { list() })
    names(pos) <- names(type.A) <- names(type.B) <- true.widths
    
	# Single chromosome, for convenience.
	distances <- round(runif(npeaks * width.n, 3000, 6000))
    
    for(width.i in seq_len(width.n)) {
    	chosen.drop.A <- sample(seq_len(6), nde, replace=TRUE)
        chosen.drop.B <- sample(seq_len(6), nde, replace=TRUE)
        
        for(section in seq_len(3)) {
            type.A[[width.i]][[section]][up.pk] <- bitwAnd(chosen.drop.A,
                c(0x1, 0x2, 0x4)[section]) > 0L
            type.B[[width.i]][[section]][down.pk] <- bitwAnd(chosen.drop.B,
                c(0x1, 0x2, 0x4)[section]) > 0L
        }
        
    	pos[[width.i]][[1]] <- cumsum(distances)[
            seq_len(npeaks) + (width.i - 1) * npeaks]
        if(!is.tf) {
        	pos[[width.i]][[2]] <- pos[[width.i]][[1]] +
                true.widths[width.i] / 2
        	pos[[width.i]][[3]] <- pos[[width.i]][[1]] + true.widths[width.i]
            sizes <- c(chr1 = max(sapply(pos, '[[', 3)) + 10000)
        } else {
            sizes <- c(chr1 = max(sapply(pos, '[[', 1)) + 10000)
        }
    	
    }
	chrs <- rep('chr1', npeaks)

	fnames <- list()
	for (lib in seq_len(length(grouping))) {
		fname <- paste0(all.fix, '_it', it, '_out_', lib, '.sam')
        message(paste(Sys.time(), 'creating', fname))
		for(width.i in seq_len(width.n)) {
            if(!is.tf) {
    			# Simulating complex histone marks.
    			if (grouping[lib] == 'A') {
     				drop <- type.A[[width.i]]
    			} else {
     				drop <- type.B[[width.i]]
    			}
                for(section in seq_len(3)) {
        			peakFile(fname, chrs = chrs[drop[[section]]],
                        pos = pos[[width.i]][[section]][drop[[section]]],
                        mu = base.mu, disp = disp[[width.i]][drop[[section]]],
        				sizes = sizes, fraglen = fraglen,
                        width = true.widths[width.i],
                        append = width.i != 1 | section != 1)
                }
            } else {
                # Simulating simple TF changes.
                cur.mu <- rep(base.mu, npeaks)
                cur.mu <- rep(base.mu, npeaks)
    			if (grouping[lib]== 'A') { 
    				cur.mu[down.pk] <- down.mu
    				cur.mu[up.pk] <- up.mu
    			} else {
    				cur.mu[up.pk] <- down.mu
    				cur.mu[down.pk] <- up.mu
    			}
    			peakFile(fname, chrs = chrs, pos = pos[[width.i]][[1]],
                    mu = cur.mu, disp = disp, sizes = sizes, fraglen = fraglen,
                    width = true.widths[width.i], tf = TRUE)
                }
            }
		}
		fnames[[lib]] <- fname
	}
	
	fnames <- unlist(fnames)
	addBackground(fnames, sizes = sizes, width = 2000 * 3, rlen = 10,
		dispersion = dispersion, prior.df = prior.df, append = TRUE)
    
    bam.files <- sapply(fnames, function(fname) {
        prefix <- sub("\\.sam$", "", basename(fname))
        bam <- paste0(prefix, '.bam')
        message(paste(Sys.time(), 'creating', bam))
        asBam(fname, prefix)
        return(bam)
    })
	if(opt$clean) unlink(fnames)

    lfile <- paste0(all.fix, '_it', it, '_log.txt')
    message(paste(Sys.time(), 'creating', lfile))
	for(width.i in seq_len(width.n)) {
        if(is.tf) {
            write.table(file = lfile,
                data.frame(chr = chrs[up.pk],
                    start = pos[[width.i]][[1]][up.pk] - radius,
                    end = pos[[width.i]][[1]][up.pk] + radius,
                    logFC = 1, truewidth = true.widths[width.i]),
    			row.names = FALSE, sep = '\t', quote = FALSE,
                append = width.i != 1)
            write.table(file = lfile,
                data.frame(chr = chrs[down.pk],
                    start = pos[[width.i]][[1]][down.pk] - radius,
                    end = pos[[width.i]][[1]][down.pk] + radius,
                    logFC = 1, truewidth = true.widths[width.i]),
    			row.names = FALSE, sep = '\t', quote = FALSE, append = TRUE,
                col.names = FALSE)
        } else {
            for(section in seq_len(3)) {
                ## Type A
                write.table(file = lfile, 
                    data.frame(chr = chrs[!type.A[[width.i]][[section]]],
                        start = pos[[width.i]][[section]][
                            !type.A[[width.i]][[section]]] - radiuses[width.i],
                        end = pos[[width.i]][[section]][
                            !type.A[[width.i]][[section]]] + radiuses[width.i],
                        logFC = -1, name = which(!type.A[[width.i]][[section]]),
                        truewidth = true.widths[width.i]),
                    row.names = FALSE, sep= '\t', quote = FALSE,
                    append = width.i != 1 | section != 1,
                    col.names = width.i == 1 & section == 1)
            
                ## Type B
                write.table(file = lfile, 
                    data.frame(chr = chrs[!type.B[[width.i]][[section]]],
                        start = pos[[width.i]][[section]][
                            !type.B[[width.i]][[section]]] - radiuses[width.i],
                        end = pos[[width.i]][[section]][
                            !type.B[[width.i]][[section]]] + radiuses[width.i],
                        logFC = 1, name = which(!type.B[[width.i]][[section]]),
                        truewidth = true.widths[width.i]),
                    row.names = FALSE, sep= '\t', quote = FALSE,
                    append = TRUE, col.names = FALSE)
            }
        }
        
	}

	#############################################################
	### Running MACS with DiffBind.
	#############################################################
	
	peakdir <- paste0(all.fix, '_peaks_it', it)
	dir.create(peakdir, showWarnings=FALSE)
	prefix <- sub('\\.bam$', '', basename(bam.files))
	gsize <- sum(as.numeric(sizes))
	
	for (peakcaller in c('MACS')) {
		if (peakcaller=='MACS') {
			pktype <- 'macs'
 		    for (x in seq_len(length(bam.files))) { 
				oprefix <- file.path(peakdir, paste0('macs_', prefix[x]))
                message(paste(Sys.time(), 'running Macs2 for', bam.files[x]))
				runMACS2(bam.files[x], oprefix, fraglen = fraglen,
                    gsize = gsize, format = 'BAM')
				
			}
		}
        
        all.peakfiles <-  dir(peakdir, pattern = 'peaks.xls', full.names = TRUE)
        
		# Running DiffBind.
        message(paste(Sys.time(), 'running DiffBind'))
		current <- dba(sampleSheet = data.frame(SampleID = prefix,
            Condition = grouping, bamReads = bam.files,
            Peaks = all.peakfiles, PeakCaller = pktype), minOverlap = 2)
		current <- dba.count(current, fragmentSize = fraglen,
            bRemoveDuplicates = FALSE)
		current <- dba.contrast(current, group1 = current$masks$A,
            group2 = current$masks$B, name1 = 'A', name2 = 'B', minMembers = 2)
		current <- dba.analyze(current, method = DBA_EDGER,
            bTagwise = TRUE) # Avoiding tagwise.dispersion when data is homoskedastic.
		out <- dba.report(current, th = 1, method = DBA_EDGER)
		
		curtab <- as.data.frame(mcols(out))
		names(curtab)[names(curtab)=='Fold'] <- 'logFC'
		for (cutoff in fdr.thresholds) {
			resultDump(out, curtab, cutoff, out = ofile)
			result <- assessChIP(ofile, lfile, tol = NA, checkfc = FALSE) # Don't bother checking fold change, as it isn't well defined for complex events.
			dump2file(paste('DiffBind +', peakcaller), cutoff, result)
		}
	}

	#############################################################
	### csaw, with its combined window methodology.
	#############################################################

	xparam <- readParam(dedup=FALSE)

	
	test.widths <- c(50, 150, 250)
	names <- c('csaw.short', 'csaw', 'csaw.long')

	for (w in seq_along(test.widths)) {
        message(paste(Sys.time(), 'running csaw with width', test.widths[w]))
		data <- windowCounts(bam.files, width = test.widths[w], ext = fraglen,
            param = xparam, filter = 20)
		binned <- windowCounts(bam.files, width = 2000, bin = TRUE,
            param = xparam)

		bin.ab <- scaledAverage(asDGEList(binned),
            scale = median(getWidths(binned)) / median(getWidths(data)))
		threshold <- median(bin.ab) + log2(2)
		keep <- aveLogCPM(asDGEList(data)) > threshold

		data <- data[keep, ]
		tabres <- analyzeQLCounts(assay(data), design, totals = data$totals)
		merged <- mergeWindows(rowRanges(data), tol=100, max.width=5000)
		tabneg <- combineTests(merged$id, tabres)

		for (cutoff in fdr.thresholds) { 
			resultDump(merged$region, tabneg, cutoff, out = ofile)
			result <- assessChIP(ofile, lfile, tol=NA, checkfc=FALSE) 
			dump2file(names[w], cutoff, result)
		}
	}
    
	#############################################################
	### derfinder with smoothing of binding signal
	#############################################################

    names(bam.files) <- gsub('.bam', '', bam.files)
    fullCov <- fullCoverage(files = bam.files, chrs = unique(chrs),
        chrlens = sizes, mc.cores = 1)
    names <- c('dbfinder.short', 'dbfinder', 'dbfinder.long')
    for(k in seq_along(test.widths)) {
        message(paste(Sys.time(), 'running derfinder with k =',
            test.widths[k] - 1))
        
        
        cuts <- seq(0.75, 2, by = 0.25)
        
            

        for(cut in cuts) {
            message(paste(Sys.time(), 'using cutoff', cut))
            regionMat <- regionMatrix(fullCov, maxClusterGap = 3000L, L = 10,
                cutoff = cut, returnBP = FALSE, smoothMean = TRUE, 
                smoothFun = bumphunter::runmedByCluster, k = test.widths[k] - 1,
                mc.cores = 1)
            
            regFile <- paste0('regionMatrix_', allfix, '_width', test.widths[k]-1, '_cut', cut, '.Rdata')
            message(paste(Sys.time(), 'saving', regFile))
            save(regionMat, file = regFile)
            
            message(paste(Sys.time(), 'running edgeR analysis'))
            keep <- width(regionMat$chr1$regions) > 10
            counts <- round(regionMat$chr1$coverageMatrix[keep, ], 0)
            tabres <- analyzeQLCounts(counts, design)
            for (cutoff in fdr.thresholds) {
                resultDump(regionMat$chr1$regions[keep], tabres, cutoff,
                    out = ofile)
                result <- assessChIP(ofile, lfile, tol = NA, checkfc = FALSE)
                dump2file(paste(names[k], cut), cutoff, result)
            }
        }
        
    }
    
    
################################################################################
    # Cleaning up.
    
    if(opt$clean) {
        unlink(bam.files)
        unlink(paste0(bam.files, '.bai'))
        unlink(lfile)
        unlink(ofile)
    }
    
}

################################################################################

## Reproducibility info
proc.time()
options(width = 120)
session_info()
