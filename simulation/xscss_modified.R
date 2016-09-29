## Functions taken from the xcss package available at http://bioinf.wehi.edu.au/csaw/
## with license LGPL (>= 2)
## Modifications summary: drop C dependencies for writing the SAM files
## because as-is xscss cannot be installed in our computing environment.

library('Rsamtools')
library('edgeR')
library('limma')
library('csaw')
library('S4Vectors')
library('IRanges')
library('GenomeInfoDb')


## New functions:

writeSAM <- function(fname, chr, y, width, rlen, type = 'hist', pos, max, alpha, beta) {
    ## Check that the file already exists
    stopifnot(file.exists(fname))
    
    ## Equivalent of i and j (1-based) from the C code
    i <- rep(seq_len(length(y)), y)
    j <- unlist(sapply(y, seq_len))
    
    ## Choose strand
    strand <- sample(c(TRUE, FALSE), length(i), replace = TRUE)
    
    ## Define current position: changes by type
    if(type == 'hist') {
        stopifnot(!missing(pos))
        stopifnot(!missing(max))
        stopifnot(!missing(alpha))
        stopifnot(!missing(beta))
        cur_pos <- rep(pos, y)
        cur_pos <- cur_pos - width / 2 + round(width * rbeta(length(i), alpha, beta))
        stopifnot(all(cur_pos <= max | cur_pos >= 1))
    } else if (type == 'bg') {
        message(paste(Sys.time(), 'calculating bg positions'))
        cur_pos <- round(runif(length(i)) * width, 0) + width * (i - 1) + 1
    }
    
    ## Write file
    df <- data.frame(
        'qname' = paste0('r.', chr, '.', j, '.', i),
        'flag' = ifelse(strand, 0, 16),
        'rname' = chr,
        'pos' = cur_pos,
        'mapq' = 200,
        'cigar' = paste0(rlen, 'M'),
        'rnext' = '*',
        'pnext' = 0,
        'tlen' = 0,
        'seq' = paste(rep('N', rlen), collapse = ''),
        'qual' = paste(rep('h', rlen), collapse = '')
    )
    message(paste(Sys.time(), 'writing to file', fname))
    write.table(file = fname, df, append = TRUE, quote = FALSE, row.names = FALSE, sep = '\t', col.names = FALSE)
}


## Functions modified:

peakFile <- function(fname, chrs, pos, mus, disp, sizes, fraglen = 100,
    width = 1000, rlen = 10, tf = FALSE, append = FALSE, alpha = 2, beta = 2) 
# Spiking in peaks more manually.
{
	if (!append) { .addHeader(fname, sizes) }

	pos <-as.integer(pos + 0.5)
	y <-as.integer(rnbinom(length(pos), mu = mus, size = 1/disp) +0.5)
	rlen <-as.integer(rlen + 0.5)
	width <-as.integer(width + 0.5)

	for (chr in unique(chrs)) {
		chosen <- chrs == chr
        out <- writeSAM(fname = fname, chr = chr, y = y[chosen], width = width,
            rlen = rlen, type = 'hist', pos = pos[chosen], max = sizes[[chr]],
            alpha = alpha, beta = beta)
	}
	return(invisible(NULL))
}

addBackground <- function(fnames, sizes, rlen=10, width=2000, back.mu=c(10, 50), prior.df=20, dispersion=0.05, append=FALSE) 
# Adding the same background regions to all files.
{
	if (!append) {
		for (fname in fnames) { .addHeader(fname, sizes) }
	}
	
	# Setting integers.
	width <- as.integer(width + 0.5)
	rlen <- as.integer(rlen + 0.5)

	# Setting up the background.
	for (chr in names(sizes)) {
		nbins <- as.integer(sizes[[chr]] / width)
		bmu <- runif(nbins, back.mu[1], back.mu[2])
		bdisp <- dispersion * prior.df / rchisq(nbins, df = prior.df)
		for (fname in fnames) {
            message(paste(Sys.time(), 'adding background to', fname))
			back.y <- as.integer(rnbinom(nbins, mu=bmu, size=1/bdisp)+0.5)
            out <- writeSAM(fname = fname, chr = chr, y = back.y,
                width = width, rlen = rlen, type = 'bg')
		}
	}
	return(invisible())
}

## Functions not modified from xcss version 1.1.2

.addHeader <- function(fname, sizes) {
	write.table(file=fname, data.frame("@SQ", paste("SN", names(sizes), sep=":"),
		paste("LN", as.integer(sizes), sep=":")), quote=FALSE, row.names=FALSE, sep="\t", col.names=FALSE)
	invisible(NULL)
}

analyzeQLCounts <- function(counts, design, totals=colSums(counts), norm.factors=1, offset=NULL, contrast=NULL, plot=FALSE, everything=FALSE)
# Runs a DE analysis on some counts using edgeR's quasi-likelihood framework (i.e. the
# QuasiSeq methods).	
#
# written by Aaron Lun
{
	# Setting up the counts.
	y<-DGEList(counts, lib.size=totals)
	y$samples$norm.factors<-norm.factors
	y$offset <- offset
	y<-estimateDisp(y, design)	
	fit<-glmQLFit(y, design, robust=TRUE)
	
	if (plot) { 
		plotBCV(y) 
		plotQLDisp(fit)
	}
	
	out<-glmQLFTest(fit, contrast=contrast)
	tab<-topTags(out, nrow(out))$table
	tab<-tab[order(as.integer(rownames(tab))),]
	if (!everything) { 
		return(tab)
	} else {
		return(list(y=y, fit=out))
	}
}

resultDump <- function(regions, tab, cutoff=0.05, out="out.tsv") 
# Dumping all the output into a specified output file.
{
	sig<-tab$FDR <= cutoff
	sig.bins<-regions[sig]
	write.table(data.frame(chr=as.character(seqnames(sig.bins)), start=start(sig.bins),
		end=end(sig.bins), tab[sig,]),
		file=out, row.names=FALSE, quote=FALSE, sep="\t")
}

assessChIP <- function(observed, known, tol=200, checkfc=TRUE) 
# We read in any tab-delimited file where the first three columns are 'chr', 'start', and 'end' of
# each predicted enriched region. We just check that the enriched sites fall within the specified
# range. The number of true positives, false positives and false negatives are reported. Both
# are reported as the number of regions rather than in terms of bases. All regions with centres
# within 'tolerance' of one another and in the same direction are merged.
#
# written by Aaron Lun
{
    obs<-read.table(observed, header=TRUE);
	if (checkfc) {
		up.o<-obs$logFC > 0
		if (is.null(up.o)) { stop("need a log-FC field in the table of observed sites") }
    }else { up.o<-rep(TRUE, nrow(obs))  }
	obranges<-GRanges(obs$chr, IRanges(obs$start, obs$end))
	if (!is.na(tol)) {
		out<-mergeWindows(obranges, sign=up.o, tol=tol)
		obranges<-out$region
		all.signs<-logical(length(obranges))
		all.signs[out$id]<-up.o
		up.o<-all.signs
	}
    
	# Reading in the known sites. Should have a 'logFC' field.
	kx<-read.table(known, header=TRUE);
	if (checkfc) {
		up.t<-kx$logFC > 0
		if (is.null(up.t)) { stop("need a log-FC field in the table of known sites") }
    } else { up.t<-rep(TRUE,nrow(kx)) }
	if (!nrow(kx)) { stop("no known sites to check against"); }
    kranges<-GRanges(kx[,1], IRanges(kx[,2], kx[,3]))
	if (is.null(kx$name)) { 
		kranges$name <- 1:length(kranges)
	} else {
		kranges$name <- kx$name
	}

	# Pulling out some states
	known.up <- kranges[up.t]
	known.down <- kranges[!up.t]
	u.olap<-findOverlaps(known.up, obranges[up.o])
	d.olap<-findOverlaps(known.down, obranges[!up.o])
	recall<-length(unique(known.up$name[queryHits(u.olap)]))+length(unique(known.down$name[queryHits(d.olap)]))
	overlapped<-length(unique(subjectHits(u.olap)))+length(unique(subjectHits(d.olap)))
	found<-length(obranges)

    return(list(overlap=overlapped, found=found, recall=recall));
}

runMACS2 <- function(file, outprefix, macs.path="macs2", threshold=NULL,
		gsize="mm", fraglen=NULL, cmd.only=FALSE, extra=NULL, format="BED") 
# This does much the same, with the new MACS2 peak caller. 
{
	cmds <- c(macs.path, "callpeak -t", file, "--gsize", gsize, "--keep-dup=all", "-f", format,
			"--outdir", dirname(outprefix), "-n", basename(outprefix), extra)
	if (!is.null(fraglen)) { cmds <- c(cmds, "--nomodel --extsize", fraglen) }
	if (!is.null(threshold)) { cmds <- c(cmds, "-p", threshold) }
	
	if (cmd.only) { return(paste(cmds, collapse=" ")) }
	if (system(paste(cmds, collapse=" "))) { stop("running MACS2 failed") }
	return(invisible(NULL))
}

