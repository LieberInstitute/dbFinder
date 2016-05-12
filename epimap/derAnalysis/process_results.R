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
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)


## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}


## Load full coverage info
rootdir <- '/dcl01/lieber/ajaffe/derRuns/derChIP/epimap/'
message(paste(Sys.time(), 'loading fullCov.Rdata'))
load(file.path(rootdir, 'CoverageInfo', 'fullCov.Rdata'))
print(object.size(fullCov), units = 'Mb')

## Load colsubset info
maindir <- file.path(rootdir, 'derAnalysis', paste0('run1-v1.5.38-',
    opt$histone))
load(file.path(maindir, 'colsubset.Rdata'))

## Subset full cov info
message(paste(Sys.time(), 'subsetting fullCov'))
fullCov <- lapply(fullCov, function(x) { x[, colsubset] })
print(object.size(fullCov), units = 'Mb')

message(paste(Sys.time(), 'saving fullCov subset'))
save(fullCov, file = file.path(rootdir, 'CoverageInfo', paste0('fullCov-',
    opt$histone, '.Rdata')))

## Load regions
message(paste(Sys.time(), 'loading fullRegions.Rdata'))
load(file.path(maindir, 'fullRegions.Rdata'))
    
## Get region coverage
message(paste(Sys.time(), 'obtaining region coverage'))
regionCoverage <- getRegionCoverage(fullCov, fullRegions)
print(object.size(regionCoverage), units = 'Mb')
rm(fullCov)

message(paste(Sys.time(), 'saving regionCov_all.Rdata'))
save(regionCoverage, file = file.path(maindir, 'regionCov_all.Rdata'))

## Construct coverage matrix
message(paste(Sys.time(), 'building coverageMatrix'))
coverageMatrix <- lapply(regionCoverage, colSums)
L <- 100 ## read length
coverageMatrix <- do.call(rbind, coverageMatrix) / L
print(object.size(coverageMatrix), units = 'Mb')

message(paste(Sys.time(), 'saving coverageMatrix'))
save(coverageMatrix, file = file.path(maindir, 'coverageMatrix.Rdata'))
rm(coverageMatrix)

## Calculate group means
message(paste(Sys.time(), 'calculating group mean coverage'))

load(file.path(maindir, 'groupInfo.Rdata'))
tIndexes <- split(seq_len(length(groupInfo)), groupInfo)
meanCoverage <- lapply(regionCoverage, function(x) {
	sapply(tIndexes, function(ii) rowMeans(x[,ii]))
})

message(paste(Sys.time(), 'saving meanCoverage.Rdata'))
save(meanCoverage, file = file.path(maindir, 'meanCoverage.Rdata'))


## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
