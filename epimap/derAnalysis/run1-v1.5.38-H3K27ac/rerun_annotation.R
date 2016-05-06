## Re-run annotation for chrs that had errors in this step due to the release
## of BioC 3.3
library('getopt')

library('GenomicRanges')
library('TxDb.Hsapiens.UCSC.hg19.knownGene')
library('bumphunter')

## Specify parameters
spec <- matrix(c(
	'chr', 'c', 1, 'character', 'Chromosome under analysis',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)


## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

message(paste(Sys.time(), 'Loading regions'))
load(file.path(opt$chr, 'regions.Rdata'))
message(paste(Sys.time(), 'analyzeChr: Annotating regions'))
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
genes <- annotateTranscripts(txdb = txdb)
annotation <- matchGenes(x = regions$regions, subject = genes)
save(annotation, file = file.path(opt$chr, 'annotation.Rdata'))

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
