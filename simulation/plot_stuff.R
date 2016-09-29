## The original version of this file is available from 
## http://bioinf.wehi.edu.au/csaw/
## and the simulation results are described in
## csaw: a Bioconductor package for differential binding analysis of ChIP-seq 
## data using sliding windows by Aaron Lun and Gordon K. Smyth
## Nucleic Acids Research 44, e45 (2015).
##
## Summary of modifications:
## 

# This plots the specificity and sensitivity values relevant to a given *.tsv file.

cutoff <- 0.05
options(width=120)

# Loading in the data.
for (prefix in c("tfx", "hist")) { 
	incoming <- paste0(prefix, "_result.tsv")
	if (prefix=="tfx") {
		methods <- c("csaw", "DBChIP + MACS", "DiffBind + HOMER", "DiffBind + MACS")
	} else {
		methods <- c("csaw", "DiffBind + HOMER", "DiffBind + MACS", "PePr")
	}

	ided <- read.table(incoming, sep="\t")
	ided<-ided[ided[,2]==cutoff,]

	out<-split(ided[,3], ided[,1])
	alpha<-bravo<-list()
	for (x in names(out)) {
		alpha[[x]]<-mean(out[[x]], na.rm=TRUE)
		bravo[[x]]<-sd(out[[x]], na.rm=TRUE)/sqrt(sum(!is.na(out[[x]])))
	}
	all.ofdrs <- unlist(alpha)
	all.ofdrse <- unlist(bravo)

	out<-split(ided[,4], ided[,1])
	alpha<-bravo<-list()
	for (x in names(out)) {
		out[[x]] <- out[[x]]/1000 # Dividing by the true number of sites.
		alpha[[x]]<-mean(out[[x]], na.rm=TRUE)
		bravo[[x]]<-sd(out[[x]], na.rm=TRUE)/sqrt(sum(!is.na(out[[x]])))
	}
	all.pows <- unlist(alpha)
	all.powse <- unlist(bravo)

	print(prefix)
	print(all.ofdrs)
	print(all.ofdrse)

	stopifnot(identical(names(all.ofdrs), names(all.pows)))
	all.names <- names(all.ofdrs)
	wo <- !grepl("\\+", all.names)
	all.names <- sub(" \\+", "\n +", all.names)
	all.names[wo] <- paste0(all.names[wo], "\n")

	for (mode in 1:2) { 
		if (mode==1L) { 
			device <- pdf
			suffix <- "pdf"
		} else {
			setEPS()
			device <- postscript
			suffix <- "eps"
		}
			
		# Plotting the observed FDR at the cutoff
		chosen <- names(all.ofdrs) %in% methods
		temp <- all.ofdrs[chosen]
		names(temp) <- " "
		temp.plus <- temp + all.ofdrse[chosen]
		
		device(paste0(prefix, "_ofdr.", suffix), width=ifelse(prefix=="tfx", 6, 4.5), height=6)
		upper <- max(temp)+0.01
		barx <- barplot(temp, col="darkgrey", ylim=c(0,upper), axis.lty=1, ylab="Observed FDR", 
			cex.lab=1.2, cex.names=1.2)
		text(barx, par("usr")[3] - upper*0.03, labels=all.names[chosen], pos=1, xpd=TRUE, cex=1)

		segments(barx, temp, barx, temp.plus, lwd=2, lend=1)
		segments(barx-0.1, temp.plus, barx+0.1, temp.plus, lwd=2)
		abline(h=cutoff, lty=3, lwd=2)
		dev.off()

		# Plotting the observed power (i.e., the number of detected regions).
		if (mode==2L) { setEPS() }

		temp <- all.pows[chosen]
		names(temp) <- " "
		temp.plus <- temp + all.powse[chosen]

		device(paste0(prefix, "_power.", suffix), width=ifelse(prefix=="tfx", 6, 4.5), height=6)
		barx <- barplot(temp, ylim=c(0, 1), axis.lty=1, ylab="Detection power", cex.lab=1.2)
		text(barx, par("usr")[3] - 0.03, labels=all.names[chosen], pos=1, xpd=TRUE, cex=1)
		
		segments(barx, temp, barx, temp.plus, lwd=2, lend=1)
		segments(barx-0.1, temp.plus, barx+0.1, temp.plus, lwd=2)
		dev.off()
	}
}

# Making an ROC plot for the TF data.

all.tf <- read.table("tf_result.tsv", sep="\t", stringsAsFactors=FALSE)
aggregator <- paste0(all.tf[,1], all.tf[,2])
all.errs <- aggregate(all.tf[,3] ~ aggregator, FUN=mean)[,2]
all.pows <- aggregate(all.tf[,4] ~ aggregator, FUN=mean)[,2]/1000
all.meth <- aggregate(all.tf[,1] ~ aggregator, FUN=function(x) { unique(x) })[,2]
all.nom <- aggregate(all.tf[,2] ~ aggregator, FUN=function(x) { unique(x) })[,2]

# Selecting only those of interest.
methods <- c("csaw", "DBChIP + MACS", "DiffBind + HOMER", "DiffBind + MACS", "PePr")
keep <- all.meth %in% methods
all.meth <- all.meth[keep]
all.errs <- all.errs[keep]
all.pows <- all.pows[keep]
all.nom <- all.nom[keep]

choices <- c(0, 1, 2, 5, 6)
filled <- c(22, 21, 24, 23, 25)
referencer <- match(all.meth, methods)

for (mode in 1:2) { 
	if (mode==1L) { 
		pdf("tf_curves.pdf")
	} else {
		setEPS()
		postscript("tf_curves.eps")
	}

	plot(0, 0, type='n', xlim=c(0, 0.2), ylim=c(0, 1), cex=1.2, cex.lab=1.2, xlab="Observed FDR", ylab="Detection power")
	for (x in unique(all.meth)) {
		cur.errs <- all.errs[all.meth==x]
		cur.pow <- all.pows[all.meth==x]
		o <- order(cur.errs)
		lines(cur.errs[o], cur.pow[o], lwd=1, lty=1, col=rgb(0,0,0))
	}

	points(all.errs, all.pows, pch=filled[referencer], cex=1.5, bg="white", col="black")
	points(all.errs[all.nom==0.05], all.pows[all.nom==0.05], pch=filled[referencer][all.nom==0.05], cex=1.8, bg="grey", col="black")
	points(all.errs[all.nom==0.05], all.pows[all.nom==0.05], pch=choices[referencer][all.nom==0.05], cex=1.8)

	abline(v=0.05, lwd=1, lty=2)

	legend("bottomright", pch=choices, cex=1.2, legend=methods, pt.cex=1.5)

	dev.off()
}