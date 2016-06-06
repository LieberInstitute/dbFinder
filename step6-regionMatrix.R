library('derfinder')
library('getopt')
library('devtools')


## Specify parameters
spec <- matrix(c(
    'maindir', 'm', 1, 'character', 'Main directory',
	'cutoff', 't', 1, 'numeric', 'Cutoff to use',
	'readLen', 'r', 1, 'integer', 'Read length',
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
maindir <- opt$maindir
cutoff <- opt$cutoff
readLen <- opt$readLen
chr <- paste0('chr', c('Y', 1:22, 'X'))[Sys.getenv('SGE_TASK_ID')]

message(Sys.time())
timeinfo <- NULL
timeinfo <- c(timeinfo, list(Sys.time()));

## Load data
load(file.path(maindir, 'CoverageInfo', paste0(chr, 'CovInfo.Rdata')))

## Load colsubset
load(file.path('..', 'derAnalysis', paste0('run1-v1.5.38-', opt$histone), 'colsubset.Rdata'))

fullCov <- list(get(paste0(chr, 'CovInfo')))
names(fullCov) <- chr
## Apply subset
fullCov[[1]]$coverage <- fullCov[[1]]$coverage[, colsubset]
timeinfo <- c(timeinfo, list(Sys.time()))
proc.time()
message(Sys.time())

## run regionMatrix
regionMat <- regionMatrix(fullCov, maxClusterGap = 3000L, L = readLen,
    cutoff = cutoff, returnBP = FALSE, smooth = TRUE, minNum = 100,
    bpSpan = 300, minInSpan = 100)
timeinfo <- c(timeinfo, list(Sys.time()))

## Save results
save(regionMat, file=paste0('regionMat-', opt$histone, '-cut', cutoff, '-', chr, '.Rdata'))
timeinfo <- c(timeinfo, list(Sys.time()))

## Save time information
save(timeinfo, file=paste0('timeinfo-', opt$histone, '-cut', cutoff, '-', chr, '.Rdata'))

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
