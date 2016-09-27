# This simulation generates a mixture of complex peak structures.
# This is done by treating each peak as a mixture of three subpeaks.
# For DB peaks, one or more of those subpeaks are randomly chosen and eliminated.

require(xscss)
require(csaw)
require(edgeR)
fraglen <- 100
iters <- 10
npeaks <- 20000
nde <- 500

prior.df <- 20
dispersion <- 1/prior.df # For the sake of possible flexibility later; right now, it just cancels out.
grouping <- c("A","A","B","B")
true.width <- 500
base.mu <- 30

if (is.tf) { 
	all.fix <- "tf"
	radius <- fraglen
	if (is.homo) { 
		prior.df <- 1e8 
		all.fix <- "tfx"
		base.mu <- 10
		up.mu <- 20
		down.mu <- 0
	} else {
		up.mu <- 45
		down.mu <- 15
	}
} else {
	is.homo <- FALSE
	all.fix <- "hist"
	radius <- true.width/2L
}

####################################################################################################

design<-model.matrix(~factor(grouping))
ofile <- paste0(all.fix, "_out.tsv")
fdr.thresholds<-c(0.01, 0.05, 0.1, 0.15, 0.2)

result.file <- paste0(all.fix, "_result.tsv")
unlink(result.file)
dump2file <- function(id, cutoff, result) {
	write.table(file=result.file, data.frame(id, cutoff, 1-result$overlap/result$found, result$recall), 
		sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
}

set.seed(21347234)

####################################################################################################

for (it in 1:iters) {
	### Generating simulated data for histone mark data.
	up.pk <- 1:nde
	down.pk <- 1:nde + nde
	disp <- prior.df*dispersion/rchisq(npeaks, df=prior.df)

	if (!is.tf) { 
		type.A.1 <- type.A.2 <- type.A.3 <- !logical(npeaks)
		chosen.drop.A <- sample(1:6, nde, replace=TRUE)
		type.A.1[up.pk] <- bitwAnd(chosen.drop.A, 0x1) > 0L
		type.A.2[up.pk] <- bitwAnd(chosen.drop.A, 0x2) > 0L
		type.A.3[up.pk] <- bitwAnd(chosen.drop.A, 0x4) > 0L

		type.B.1 <- type.B.2 <- type.B.3 <- !logical(npeaks)
		chosen.drop.B <- sample(1:6, nde, replace=TRUE)
		type.B.1[down.pk] <- bitwAnd(chosen.drop.B, 0x1) > 0L
		type.B.2[down.pk] <- bitwAnd(chosen.drop.B, 0x2) > 0L
		type.B.3[down.pk] <- bitwAnd(chosen.drop.B, 0x4) > 0L
	}

	# Single chromosome, for convenience.
	distances<-round(runif(npeaks, 10000, 20000))
	pos.1 <- cumsum(distances)	
	if (!is.tf) { 
		pos.2 <- pos.1 + true.width/2
		pos.3 <- pos.1 + true.width
		sizes <- c(chrA=max(pos.3) + 10000)
	} else {
		sizes <- c(chrA=max(pos.1) + 10000)
	}
	chrs <- rep("chrA", npeaks)

	fnames<-list()
	for (lib in 1:length(grouping)) {
		fname <- paste0(all.fix, "_out_", lib, ".sam")
		if (!is.tf) {
			# Simulating complex histone marks.
			if (grouping[lib]=="A") {
 				drop.1 <- type.A.1
				drop.2 <- type.A.2
				drop.3 <- type.A.3
			} else {
 				drop.1 <- type.B.1
				drop.2 <- type.B.2
				drop.3 <- type.B.3
			}
			peakFile(fname, chrs=chrs[drop.1], pos=pos.1[drop.1], mu=base.mu, disp=disp[drop.1],
				sizes=sizes, fraglen=fraglen, width=true.width, tf=FALSE)
			peakFile(fname, chrs=chrs[drop.2], pos=pos.2[drop.2], mu=base.mu, disp=disp[drop.2],
				sizes=sizes, fraglen=fraglen, width=true.width, tf=FALSE, append=TRUE)
			peakFile(fname, chrs=chrs[drop.3], pos=pos.3[drop.3], mu=base.mu, disp=disp[drop.3],
				sizes=sizes, fraglen=fraglen, width=true.width, tf=FALSE, append=TRUE)

		} else {
			# Simulating simple TF changes.
			cur.mu <- rep(base.mu, npeaks)
			if (grouping[lib]=="A") { 
				cur.mu[down.pk] <- down.mu
				cur.mu[up.pk] <- up.mu
			} else {
				cur.mu[up.pk] <- down.mu
				cur.mu[down.pk] <- up.mu
			}
			peakFile(fname, chrs=chrs, pos=pos.1, mu=cur.mu, disp=disp,
				sizes=sizes, fraglen=fraglen, width=true.width, tf=TRUE)
		}
		fnames[[lib]]<-fname
	}
	
	fnames<-unlist(fnames)
	addBackground(fnames, sizes=sizes, width=2000, rlen=10,
		dispersion=dispersion, prior.df=prior.df, append=TRUE)
	bam.files<-crunch2BAM(fnames)
	unlink(fnames)

    lfile <- paste0(all.fix, "_log.txt")
	if (is.tf) { 
		write.table(file=lfile, data.frame(chr=chrs[up.pk], start=pos.1[up.pk]-radius, end=pos.1[up.pk]+radius, logFC=1),
			row.names=FALSE, sep="\t", quote=FALSE)
		write.table(file=lfile, data.frame(chrs[down.pk], pos.1[down.pk]-radius, pos.1[down.pk]+radius, logFC=-1),
			row.names=FALSE, sep="\t",  quote=FALSE, append=TRUE, col.names=FALSE)
	} else {
		write.table(file=lfile, data.frame(chr=chrs[!type.A.1], start=pos.1[!type.A.1]-radius, end=pos.1[!type.A.1]+radius, logFC=-1, 
			name=which(!type.A.1)), row.names=FALSE, sep="\t", quote=FALSE)
		write.table(file=lfile, data.frame(chr=chrs[!type.A.2], start=pos.2[!type.A.2]-radius, end=pos.2[!type.A.2]+radius, logFC=-1, 
			name=which(!type.A.2)), row.names=FALSE, sep="\t", quote=FALSE, append=TRUE, col.names=FALSE)
		write.table(file=lfile, data.frame(chr=chrs[!type.A.3], start=pos.3[!type.A.3]-radius, end=pos.3[!type.A.3]+radius, logFC=-1, 
			name=which(!type.A.3)), row.names=FALSE, sep="\t", quote=FALSE, append=TRUE, col.names=FALSE)
		write.table(file=lfile, data.frame(chr=chrs[!type.B.1], start=pos.1[!type.B.1]-radius, end=pos.1[!type.B.1]+radius, logFC=1, 
			name=which(!type.B.1)), row.names=FALSE, sep="\t", quote=FALSE, append=TRUE, col.names=FALSE)
		write.table(file=lfile, data.frame(chr=chrs[!type.B.2], start=pos.2[!type.B.2]-radius, end=pos.2[!type.B.2]+radius, logFC=1, 
			name=which(!type.B.2)), row.names=FALSE, sep="\t", quote=FALSE, append=TRUE, col.names=FALSE)
		write.table(file=lfile, data.frame(chr=chrs[!type.B.3], start=pos.3[!type.B.3]-radius, end=pos.3[!type.B.3]+radius, logFC=1, 
			name=which(!type.B.3)), row.names=FALSE, sep="\t", quote=FALSE, append=TRUE, col.names=FALSE)
	}

	#############################################################
	### Running HOMER or MACS with DiffBind.
	#############################################################
	
	peakdir <- paste0(all.fix, "_peaks")
	dir.create(peakdir, showWarnings=FALSE)
	prefix <- sub("\\.bam$", "", basename(bam.files))
	gsize <- sum(as.numeric(sizes))
	
	for (peakcaller in c("HOMER", "MACS", "SICER")) {
		all.peakfiles <- list()

		if (peakcaller=="HOMER") {
			pktype <- "bed"
			hpdir <- file.path(peakdir, "homerpool")
			dir.create(hpdir)
			for (x in 1:length(bam.files)) { 
				pooldir <- file.path(hpdir, prefix[x]) 
				dir.create(pooldir)
				system(paste("makeTagDirectory", pooldir, bam.files[x], "-format sam -keepAll"))

				outers <- file.path(peakdir, paste0("homer_", prefix[x], ".txt"))
				system(paste("findPeaks", pooldir, "-style", ifelse(is.tf, "factor", "histone"), "-o", outers, "-gsize", gsize, "-fragLength", fraglen, "-tbp", 0))

				# Converting to proper BED format.
				blah <- read.table(outers)
				write.table(file=outers, blah[,c(2,3,4,1,8)], row.names=FALSE, quote=FALSE, sep="\t", col.names=FALSE)
				all.peakfiles[[x]] <- outers
			}
			unlink(hpdir, recursive=TRUE)

		} else if (peakcaller=="MACS") {
			pktype <- "macs"
 		    for (x in 1:length(bam.files)) { 
				oprefix <- file.path(peakdir, paste0("macs_", prefix[x]))
				runMACS2(bam.files[x], oprefix, fraglen=fraglen, gsize=gsize, format="BAM")
				all.peakfiles[[x]] <- paste0(oprefix, "_peaks.xls")
			}
			macs.peakfiles <- all.peakfiles # For use with DBChIP.

		} else if (peakcaller=="SICER") {
			if (is.tf || it > 1L) { next } # Only doing it for histone mark data, and just once, for demonstration.
			pktype <- "bed"
			scdir <- file.path(peakdir, "sicerbed")
			dir.create(scdir)
			if (!file.exists("~/tmp")) { dir.create("~/tmp") } 
			win <- 200
			gap <- win*2L
			eval <- 1

			all.peakfiles <- list()
			for (x in 1:length(bam.files)) {
				cur.bed <- "out.bed"
 	 	   		convertBamToBed(bam.files[x], file.path(scdir, cur.bed))
				system(paste("methods/SICER_V1.1/SICER/SICER-rb.sh", scdir, cur.bed, scdir, "Simulated", 
						2147483647, win, fraglen, 1, gap, eval))
				newname <- file.path(peakdir, paste0("sicer_", x, ".bed"))
				blah <- read.table(file.path(scdir, sprintf("out-W%i-G%i-E%i.scoreisland", win, gap, eval)))
				write.table(data.frame(blah[,1:3], paste0("peak-", 1:nrow(blah)), blah[,4]), file=newname, 
					sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
				all.peakfiles[[x]] <- newname
			}

			unlink("chr.list")
			unlink(scdir, recursive=TRUE)
		}

		# Running DiffBind.
		require(DiffBind)
		current <- dba(sampleSheet=data.frame(SampleID=prefix, Condition=grouping, 
			bamReads=bam.files, Peaks=unlist(all.peakfiles), PeakCaller=pktype), minOverlap=2)
		current <- dba.count(current, fragmentSize=fraglen, bRemoveDuplicates=FALSE)
		current <- dba.contrast(current, group1=current$masks$A, group2=current$masks$B, name1="A", name2="B", minMembers=2)
		current <- dba.analyze(current, method=DBA_EDGER, bTagwise=!is.homo) # Avoiding tagwise.dispersion when data is homoskedastic.
		out <- dba.report(current, th=1)
		
		curtab <- as.data.frame(elementMetadata(out))
		names(curtab)[names(curtab)=="Fold"] <- "logFC"
		for (cutoff in fdr.thresholds) { 
			resultDump(out, curtab, cutoff, out=ofile)
			result<-assessChIP(ofile, lfile, tol=NA, checkfc=FALSE) # Don't bother checking fold change, as it isn't well defined for complex events.
			dump2file(paste("DiffBind +", peakcaller), cutoff, result)
		}
	}

	#############################################################
	# Running MACS with DBChIP.
	#############################################################

	if (is.tf) {
		all.bed <- list()
		for (x in 1:length(bam.files)) {
			curbed <- file.path(peakdir, sub("\\.bam$", ".bed", basename(bam.files[x])))
			convertBamToBed(bam.files[x], curbed)
			all.bed[[x]] <- curbed		   
		}
		all.bed <- unlist(all.bed)
		conds <- factor(as.character(grouping))
		
#		# Re-running MACS on pooled reads for each condition.
#		pool.1 <- file.path(peakdir, "pooledA.bed")
#		system(paste(c("cat", all.bed[grouping=="A"], ">", pool.1), collapse=" "))
#		opref.1 <- file.path(peakdir, "macs_A")
#		runMACS2(pool.1, opref.1, fraglen=fraglen, gsize=gsize)
#
#		pool.2 <- file.path(peakdir, "pooledB.bed")
#		system(paste(c("cat", all.bed[grouping=="B"], ">", pool.2), collapse=" "))
#		opref.2 <- file.path(peakdir, "macs_B")
#		runMACS2(pool.2, opref.2, fraglen=fraglen, gsize=gsize)
#
#		binding.site.list <- list()
#		binding.site.list[["A"]] <- read.table(paste0(opref.1, "_peaks.xls"), header=TRUE)[,c("chr", "abs_summit", "X.log10.pvalue.")]
#		binding.site.list[["B"]] <- read.table(paste0(opref.2, "_peaks.xls"), header=TRUE)[,c("chr", "abs_summit", "X.log10.pvalue.")]
#		colnames(binding.site.list$A) <- colnames(binding.site.list$B) <- c("chr", "pos", "weight")

		# Loading the peaks, using the first peak set for each condition.
		binding.site.list <- list()
		for (x in unique(grouping)) { 
			binding.site.list[[x]] <- read.table(macs.peakfiles[[which(grouping==x)[1]]], header=TRUE)[,c("chr", "abs_summit", "X.log10.pvalue.")]
			colnames(binding.site.list[[x]]) <- c("chr", "pos", "weight")
		}

		# Running DBCHIP.
		require(DBChIP)
		bs.list <- read.binding.site.list(binding.site.list)
		consensus.site <- site.merge(bs.list)
		chip.data.list <- as.list(all.bed)
		names(chip.data.list) <- names(conds) <- paste0("x", 1:length(conds)) # Just to avoid confusion.
		dat <- load.data(chip.data.list=chip.data.list, conds=conds, consensus.site=consensus.site, data.type="BED", 
			chr.vec=names(sizes), chr.len.vec=sizes, frag.len=fraglen)
		dat <- get.site.count(dat)
		dat <- test.diff.binding(dat) # Default uses common dispersion, exactTest; no need for special behaviour when is.homo=TRUE.
		rept <- report.peak(dat, n=Inf)

		out <- GRanges(rept$chr, IRanges(rept$pos, rept$pos))
		curtab <- data.frame(logFC=log2(rept$FC.B), PValue=rept$pval, FDR=rept$FDR)
		for (cutoff in fdr.thresholds) { 
			resultDump(out, curtab, cutoff, out=ofile)
			result<-assessChIP(ofile, lfile, tol=NA, checkfc=FALSE)
			dump2file("DBChIP + MACS", cutoff, result)
		}

		unlink(all.bed)
	}

	#############################################################
	# Running PePr.
	#############################################################

	outname <- file.path(peakdir, "pepr") 
	first.set <- paste(bam.files[grouping=="A"], collapse=",")
	second.set <- paste(bam.files[grouping=="B"], collapse=",")
	system(sprintf("python methods/PePr-master/PePr/PePr.py -c %s --chip2 %s -n %s -f BAM --peaktype=%s --diff -s %i", 
		first.set, second.set, outname, ifelse(is.tf, "sharp", "broad"), fraglen/2L))

	collected.up <- read.table(paste0(outname, "__PePr_up_peaks.bed"))
	collected.down <- read.table(paste0(outname, "__PePr_down_peaks.bed"))
	all.collected <- c(GRanges(collected.up[,1], IRanges(collected.up[,2], collected.up[,3])),
			GRanges(collected.down[,1], IRanges(collected.down[,2], collected.down[,3])))
	curtab <- rbind(collected.up[,5:6], collected.down[,5:6])
	curtab <- cbind(c(rep(1, nrow(collected.up)), rep(-1, nrow(collected.down))), curtab)
	colnames(curtab) <- c("logFC", "PValue", "FDR")
	
	for (cutoff in fdr.thresholds) { 
		# We can't turn up the default window threshold (otherwise everything gets aggregated into a peak).
		# We'll have to work with the FDRs that are reported, most of which are below 0.05. As such, 
		# results will be conservative for most of the cut-offs.
		resultDump(all.collected, curtab, cutoff, out=ofile)
		result<-assessChIP(ofile, lfile, tol=NA, checkfc=FALSE)
		dump2file("PePr", cutoff, result)
	}

	unlink(list.files(pattern="debug\\.log$"))

	#############################################################
	# Running multiGPS, but only once for each data set.
	#############################################################
	
	if (is.tf && it==1L) { 
		gpsdir <- file.path(peakdir, "multigps") 
		dir.create(gpsdir)
		file.spec <- paste0('--expt', grouping, "-", as.integer(duplicated(grouping))+1, " ", bam.files)
		gdata <- file.path(gpsdir, "gdata.data")
		write.table(data.frame("chrA", gsize), file=gdata, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

		system(sprintf("java -Xmx4G -jar methods/multigps_v0.5.jar --geninfo %s %s --out %s", 
			gdata, paste(file.spec, collapse=" "), gpsdir))

		file.copy(file.path(gpsdir, "multigps_A.events"), file.path(peakdir, "multigps_A.events"))
		file.copy(file.path(gpsdir, "multigps_B.events"), file.path(peakdir, "multigps_B.events"))
		null <- data.frame(character(0), matrix(0, ncol=8, nrow=0))
		collected.up <- tryCatch({
			read.table(file.path(peakdir, "multigps_A.events"), stringsAsFactors=FALSE)
		}, error=function(w) {
			null
		})
		collected.down <- tryCatch({
			read.table(file.path(peakdir, "multigps_B.events"), stringsAsFactors=FALSE)
		}, error=function(w) { 
			null
		})

		cur.str <- c(collected.up[,1], collected.down[,1])
		cur.pos <- as.integer(sub(".*:", "", cur.str))
		cur.chr <- sub(":.*", "", cur.str)
		granges <- GRanges(cur.chr, IRanges(cur.pos, cur.pos))
		curtab <- data.frame(logCPM=c(collected.up[,6], collected.down[,6]), 
			logFC=c(collected.up[,7], -collected.down[,7]),
			FDR=2^c(collected.up[,8], collected.down[,8])) 
			# Assuming already BH-adjusted (see setting of DEPVal in EdgeRDifferentialEnrichment.java,
			# ultimately used in setInterCondP to set interCondP for final printing to *.events via
			# writeBindingEventFiles in DifferentialTester.java.

		for (cutoff in fdr.thresholds) { 
			resultDump(granges, curtab, cutoff, out=ofile)
			result<-assessChIP(ofile, lfile, tol=NA, checkfc=FALSE) 
			dump2file("multiGPS", cutoff, result)
		}

		unlink(gdata)
		unlink(gpsdir, recursive=TRUE)
	}

#	regs <- read.table(outers, sep="\t")
#	regions <- GRanges(regs[,2], IRanges(regs[,3], regs[,4]))
#	data<-regionCounts(bam.files, regions, ext=1L, param=readParam(dedup=FALSE))
#	keep<-rowSums(assay(data))>=filter
#	counts<-assay(data)[keep,,drop=FALSE]
#	tabhom<-analyzeQLCounts(counts, design, totals=data$totals)
#	unlink(pooldir, recursive=TRUE)

	#############################################################
	### csaw, with its combined window methodology.
	############################################################

	xparam <- readParam(dedup=FALSE)

	if (is.tf) { 
		test.widths <- 10
		names <- c("csaw")
	} else {
		test.widths <- c(50, 150, 250)
		names <- c("csaw.short", "csaw", "csaw.long")
	}

	for (w in seq_along(test.widths)) { 
		data <- windowCounts(bam.files, width=test.widths[w], ext=fraglen, param=xparam, filter=20)
		binned <- windowCounts(bam.files, width=2000, bin=TRUE, param=xparam)

		bin.ab <- scaledAverage(asDGEList(binned), scale=median(getWidths(binned))/median(getWidths(data)))
		threshold <- median(bin.ab) + log2(2)
		keep <- aveLogCPM(asDGEList(data)) >  threshold

		data <- data[keep,]
		tabres <- analyzeQLCounts(assay(data), design, totals=data$totals)
		merged <- mergeWindows(rowRanges(data), tol=100, max.width=5000)
		tabneg <- combineTests(merged$id, tabres)

		for (cutoff in fdr.thresholds) { 
			resultDump(merged$region, tabneg, cutoff, out=ofile)
			result<-assessChIP(ofile, lfile, tol=NA, checkfc=FALSE) 
			dump2file(names[w], cutoff, result)
		}
	}
}

####################################################################################################
# Cleaning up.

unlink(bam.files)
unlink(paste0(bam.files, ".bai"))
unlink(lfile)
unlink(ofile)

####################################################################################################
