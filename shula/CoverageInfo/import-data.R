## Generate fullCov object from BED files

library('rtracklayer')
library('BiocParallel')
library('devtools')
library('BSgenome.Hsapiens.UCSC.hg19')

## Define files
beds <- dir('/dcs01/ajaffe/ChIPseq/Shulha2013/BED', pattern = 'c', full.names = TRUE)
names(beds) <- dir('/dcs01/ajaffe/ChIPseq/Shulha2013/BED', pattern = 'c')

## Load data for all chrs per file
system.time(
    bedGR <- bplapply(beds, function(bed) {
        library('rtracklayer')
        b <- BEDFile(bed)
        import(b)
    }, BPPARAM = SnowParam(workers = 10))
)
print(object.size(bedGR), units = 'Gb')

## Set the chr lengths
bedGR <- lapply(bedGR, function(gr) {
    seqs <- names(seqlengths(gr))
    match.names <- match(seqs, names(seqlengths(BSgenome.Hsapiens.UCSC.hg19)))
    match.names <- match.names[!is.na(match.names)]
    seqlengths(gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[match.names]
    return(gr)
})

## Calculate coverage
system.time(
    bedCov <- bplapply(bedGR, function(bed) {
        library('GenomicRanges')
        coverage(bed)
    }, BPPARAM = SnowParam(workers = 10))
)
print(object.size(bedCov), units = 'Gb')

## Build fullCov object
chrs <- paste0('chr', c(1:22, 'X', 'Y', 'M'))
system.time(
    fullCov <- lapply(chrs, function(chr) {
        DataFrame(lapply(bedCov, '[[', chr))
    })
)
names(fullCov) <- chrs

## Save the coverage data
save(fullCov, file = 'fullCov.Rdata')

## Filter the data and save it by chr
myFilt <- function(chr, rawData, cutoff) {
    library('derfinder')
    message(paste(Sys.time(), 'Filtering chromosome', chr))
    
	## Filter the data
	res <- filterData(data = rawData, cutoff = cutoff, index = NULL)
	
	## Save it in a unified name format
	varname <- paste0(chr, 'CovInfo')
	assign(varname, res)
	output <- paste0(varname, '.Rdata')
	
	## Save the filtered data
	save(list = varname, file = output, compress='gzip')
	
	## Finish
	return(invisible(NULL))
}

message(paste(Sys.time(), 'Filtering and saving the data with cutoff', 2))
filteredCov <- bpmapply(myFilt, names(fullCov), fullCov, BPPARAM = SnowParam(workers = 10), MoreArgs = list(cutoff = 2))

source('check-filter.R')

## Reproducibility
proc.time()
options(width = 120)
session_info()
