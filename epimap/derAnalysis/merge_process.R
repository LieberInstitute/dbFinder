## Merge chr results for H3K27ac

library('derfinder')
library('GenomicRanges')
library('devtools')

rootdir <- '/dcl01/lieber/ajaffe/derRuns/derChIP/epimap/'
maindir <- file.path(rootdir, 'derAnalysis', 'run1-v1.5.38-H3K27ac')

## Merge coverage matrix
chrs <- paste0('chr', c(1:22, 'X', 'Y')) # chrM has no regions

loadCov <- function(chr) {
    load(file.path(maindir, paste0('coverageMatrix', chr, '.Rdata')))
    return(coverageMatrix)
}

message(paste(Sys.time(), 'building coverageMatrix'))
coverageMatrix <- lapply(chrs, loadCov)
coverageMatrix <- do.call(rbind, coverageMatrix)
print(object.size(coverageMatrix), units = 'Mb')

message(paste(Sys.time(), 'saving coverageMatrix'))
save(coverageMatrix, file = file.path(maindir, 'coverageMatrix.Rdata'))
rm(coverageMatrix)

loadMean <- function(chr) {
    load(file = file.path(maindir, paste0('meanCoverage-', chr,
        '.Rdata')))
    return(meanCoverage)
}

message(paste(Sys.time(), 'building mean coverage'))
meanCoverage <- lapply(chrs, loadMean)
if(class(meanCoverage[[1]]) == 'list') meanCoverage <- unlist(meanCoverage)
print(object.size(meanCoverage), units = 'Mb')

message(paste(Sys.time(), 'saving meanCoverage.Rdata'))
save(meanCoverage, file = file.path(maindir, 'meanCoverage.Rdata'))
rm(meanCoverage)

loadReg <- function(cr) {
    load(file = file.path(maindir, paste0('regionCov-', chr,
        '.Rdata')))
    return(regionCoverage)
}

message(paste(Sys.time(), 'obtaining region coverage'))
regionCoverage <- lapply(chrs, loadReg)
if(class(regionCoverage[[1]]) == 'list') {
    regionCoverage <- unlist(regionCoverage)
}
print(object.size(regionCoverage), units = 'Mb')

message(paste(Sys.time(), 'saving regionCov_all.Rdata'))
save(regionCoverage, file = file.path(maindir, 'regionCov_all.Rdata'))

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
