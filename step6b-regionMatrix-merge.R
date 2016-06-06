library('derfinder')
library('GenomicRanges')
library('getopt')
library('devtools')


## Specify parameters
spec <- matrix(c(
	'cutoff', 't', 1, 'numeric', 'Cutoff to use',
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


## Input options used
cutoff <- opt$cutoff

## Load data
chrs <- paste0('chr', c(1:22, 'X', 'Y'))
regionMat <- lapply(chrs, function(chr) {
    message(paste(Sys.time(), 'processing', chr))
    load(paste0('regionMat-', opt$histone, '-cut', cutoff, '-', chr, '.Rdata'))
    res <- regionMat
    return(res)
})

## Merge
regionMat <- do.call(c, regionMat)

## Save
message(paste(Sys.time(), 'saving full results'))
save(regionMat, file = paste0('regionMat-', opt$histone, '-cut', cutoff, '.Rdata'))

## Extract info
message(paste(Sys.time(), 'extracting regions'))
fullRegions  <- unlist(GRangesList(lapply(regionMat, '[[', 'regions')))
names(fullRegions) = NULL

message(paste(Sys.time(), 'saving regions'))
save(fullRegions, file = paste0('fullRegions-', opt$histone, '-cut', cutoff, '.Rdata'))

message(paste(Sys.time(), 'extracting coverage matrix'))
coverageMatrix <- do.call("rbind", lapply(regionMat, '[[', 'coverageMatrix'))

message(paste(Sys.time(), 'saving coverage matrix'))
save(coverageMatrix, file = paste0('coverageMatrix-', opt$histone, '-cut', cutoff, '.Rdata'))

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
