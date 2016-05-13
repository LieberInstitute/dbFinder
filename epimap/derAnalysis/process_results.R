## Process derfinder results

## Load libraries
library('getopt')

## Available at http://www.bioconductor.org/packages/release/bioc/html/derfinder.html
library('derfinder')
library('GenomicRanges')
library('devtools')

## Specify parameters
spec <- matrix(c(
    'histone', 'i', 1, 'character', 'For epimap, the histone mark to use. Either H3K27ac or H3K4me3.',
    'chr', 'c', 2, 'character', 'Chromosome number',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)


## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

chrFlag <- !is.null(opt$chr) | opt$chr != 'chrall'

## Load colsubset info
rootdir <- '/dcl01/lieber/ajaffe/derRuns/derChIP/epimap/'
maindir <- file.path(rootdir, 'derAnalysis', paste0('run1-v1.5.38-',
    opt$histone))
load(file.path(maindir, 'colsubset.Rdata'))

## Load regions
message(paste(Sys.time(), 'loading fullRegions.Rdata'))
load(file.path(maindir, 'fullRegions.Rdata'))

if(!chrFlag) {
    ## Load full coverage info
    message(paste(Sys.time(), 'loading fullCov.Rdata'))
    load(file.path(rootdir, 'CoverageInfo', 'fullCov.Rdata'))
    print(object.size(fullCov), units = 'Mb')    

    ## Subset full cov info
    message(paste(Sys.time(), 'subsetting fullCov'))
    fullCov <- lapply(fullCov, function(x) { x[, colsubset] })
    print(object.size(fullCov), units = 'Mb')

    message(paste(Sys.time(), 'saving fullCov subset'))
    save(fullCov, file = file.path(rootdir, 'CoverageInfo', paste0('fullCov-',
        opt$histone, '.Rdata')))
} else {
    chr <- paste0('chr', opt$chr)
    
    ## Load coverage info
    message(paste(Sys.time(), 'loading fullCov.Rdata'))
    load(file.path(rootdir, 'CoverageInfo', paste0(chr, 'CovInfo.Rdata')))
    eval(parse(text=paste0('fullCov <- list(', chr, ' = ', chr, 'CovInfo$coverage)')))
    eval(parse(text=paste0('rm(', chr, 'CovInfo)')))
    print(object.size(fullCov), units = 'Mb')

    ## Subset full cov info
    message(paste(Sys.time(), 'subsetting fullCov'))
    fullCov <- lapply(fullCov, function(x) { x[, colsubset] })
    print(object.size(fullCov), units = 'Mb')
    
    message(paste(Sys.time(), 'saving fullCov subset'))
    save(fullCov, file = file.path(rootdir, 'CoverageInfo', paste0('fullCov-',
        chr, '-', opt$histone, '.Rdata')))
        
    ## Subset regions
    fullRegions <- fullRegions[seqnames(fullRegions) == chr]
}
    
## Get region coverage
message(paste(Sys.time(), 'obtaining region coverage'))
regionCoverage <- getRegionCoverage(fullCov, fullRegions)
print(object.size(regionCoverage), units = 'Mb')
rm(fullCov)

message(paste(Sys.time(), 'saving regionCov_all.Rdata'))
if(chrFlag) {
    save(regionCoverage, file = file.path(maindir, paste0('regionCov-', chr,
        '.Rdata')))
} else {
    save(regionCoverage, file = file.path(maindir, 'regionCov_all.Rdata'))
}


## Construct coverage matrix
message(paste(Sys.time(), 'building coverageMatrix'))
coverageMatrix <- lapply(regionCoverage, colSums)
L <- 100 ## read length
coverageMatrix <- do.call(rbind, coverageMatrix) / L
print(object.size(coverageMatrix), units = 'Mb')

message(paste(Sys.time(), 'saving coverageMatrix'))
if(chrFlag) {
    save(coverageMatrix, file = file.path(maindir, paste0('coverageMatrix',
        chr, '.Rdata')))
} else {
    save(coverageMatrix, file = file.path(maindir, 'coverageMatrix.Rdata'))
}

rm(coverageMatrix)

## Calculate group means
message(paste(Sys.time(), 'calculating group mean coverage'))

load(file.path(maindir, 'groupInfo.Rdata'))
tIndexes <- split(seq_len(length(groupInfo)), groupInfo)
meanCoverage <- lapply(regionCoverage, function(x) {
	sapply(tIndexes, function(ii) rowMeans(x[,ii]))
})
print(object.size(meanCoverage), units = 'Mb')

message(paste(Sys.time(), 'saving meanCoverage.Rdata'))
if(chrFlag) {
    save(meanCoverage, file = file.path(maindir, paste0('meanCoverage-', chr,
        '.Rdata')))
} else {
    save(meanCoverage, file = file.path(maindir, 'meanCoverage.Rdata'))
}

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
