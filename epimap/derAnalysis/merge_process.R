## Merge chr results for H3K27ac

library('derfinder')
library('GenomicRanges')
library('devtools')

rootdir <- '/dcl01/lieber/ajaffe/derRuns/derChIP/epimap/'
maindir <- file.path(rootdir, 'derAnalysis', 'run1-v1.5.38-H3K27ac')

## Merge coverage matrix
chrs <- paste0('chr', c(1:22, 'X', 'Y')) # chrM has no regions

loadCov <- function(chr) {
    load(file.path(maindir, paste0('coverageMatrix-', chr, '.Rdata')))
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

manualUnlist <- function(covList, names) {
    covList <- c(covList[[1]], covList[[2]], covList[[3]], covList[[4]],
        covList[[5]], covList[[6]], covList[[7]], covList[[8]], covList[[9]],
        covList[[10]], covList[[11]], covList[[12]], covList[[13]],
        covList[[14]], covList[[15]], covList[[16]], covList[[17]],
        covList[[18]], covList[[19]], covList[[20]], covList[[21]],
        covList[[22]], covList[[23]], covList[[24]]
    )
    names(covList) <- names
    return(covList)
}

message(paste(Sys.time(), 'loading fullRegions'))
load(file.path(maindir, 'fullRegions.Rdata'))
name_groups <- split(seq_len(length(fullRegions)), seqnames(fullRegions))
stopifnot(identical(names(name_groups), chrs))
regNames <- unlist(name_groups)

message(paste(Sys.time(), 'building mean coverage'))
meanCoverage <- lapply(chrs, loadMean)
meanCoverage <- manualUnlist(meanCoverage, regNames)
print(object.size(meanCoverage), units = 'Mb')

message(paste(Sys.time(), 'saving meanCoverage.Rdata'))
save(meanCoverage, file = file.path(maindir, 'meanCoverage.Rdata'))
rm(meanCoverage)

loadReg <- function(chr) {
    load(file = file.path(maindir, paste0('regionCov-', chr,
        '.Rdata')))
    return(regionCoverage)
}

message(paste(Sys.time(), 'obtaining region coverage'))
regionCoverage <- lapply(chrs, loadReg)
regionCoverage <- manualUnlist(regionCoverage, regNames)
print(object.size(regionCoverage), units = 'Mb')

message(paste(Sys.time(), 'saving regionCov_all.Rdata'))
save(regionCoverage, file = file.path(maindir, 'regionCov_all.Rdata'))

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
