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


## Save unfiltered chrs individually
xx <- lapply(chrs, function(chr) {
    varname <- paste0(chr, 'CovInfo')
    res <- list('coverage' = fullCov[[chr]], 'position' = NULL)
    assign(varname, res)
    output <- paste0(varname, '.Rdata')
    
    ## Save the raw data
    save(list = varname, file = output, compress = 'gzip')
    
    ## Finish
    return(invisible(NULL))
})

## Proceed to filter the coverage data
rm(fullCov)

fullCov_files <- as.list(dir(pattern = 'chr'))
names(fullCov_files) <- sapply(fullCov_files, function(x) gsub('CovInfo.Rdata', '', x))

## Filter the data and save it by chr
myFilt <- function(chr, rawData_file, cutoff, totalMapped = NULL, targetSize = 80e6) {
    library('derfinder')
    
    ## Load raw data
    message(paste(Sys.time(), 'Loading raw file', rawData_file, 'for', chr))
    load(rawData_file)
    rawData <- get(paste0(chr, 'CovInfo'))
    
	## Filter the data
    message(paste(Sys.time(), 'Filtering chromosome', chr))
	res <- filterData(data = rawData$coverage, cutoff = cutoff, index = NULL,
        totalMapped = totalMapped, targetSize = targetSize)
	
	## Save it in a unified name format
	varname <- paste0(chr, 'CovInfo')
	assign(varname, res)
	output <- paste0(varname, '-filtered.Rdata')
	
	## Save the filtered data
	save(list = varname, file = output, compress='gzip')
	
	## Finish
	return(invisible(NULL))
}


message(paste(Sys.time(), 'Filtering and saving the data with cutoff', opt$cutoff))
filteredCov <- bpmapply(myFilt, names(fullCov_files), fullCov_files, BPPARAM = SnowParam(workers = opt$mcores), MoreArgs = list(cutoff = opt$cutoff))

## Done!
proc.time()
options(width = 120)
session_info()
