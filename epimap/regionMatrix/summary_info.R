## Compute basic summary information from the results
# module load R/3.3.x
# mkdir -p logs
# Rscript summary_info.R > logs/summary_info.txt 2>&1

library('GenomicRanges')

histones <- c('H3K4me3', 'H3K27ac')
cutoffs <- c(5, 10)

regions <- lapply(histones, function(histone) {
    res <- lapply(cutoffs, function(cutoff) {
        load(file.path(paste0(histone, '-cut', cutoff), 'fullRegions.Rdata'))
        return(fullRegions)
    })
    names(res) <- cutoffs
    return(res)
})
names(regions) <- histones

## Number of ERs
print('Number of ERs')
sapply(regions, function(x) { sapply(x, function(y) { length(y) })})

## Number of DERs
print('Number of DERs')
sapply(regions, function(x) { sapply(x, function(y) { sum(y$significantFWER == 'TRUE') })})

## Percent DE
print('Percent DE')
sapply(regions, function(x) { sapply(x, function(y) { round(sum(y$significantFWER == 'TRUE')/ length(y) * 100, 2) })})

## Overlap between sets
print('Overlap between different cutoffs for each histone mark')
lapply(regions, function(x) {
    sig <- lapply(x, function(y) { y[y$significantFWER == 'TRUE' ]})
    nonsig <- lapply(x, function(y) { y[y$significantFWER != 'TRUE' ]})
    counts <- list(
        '5 versus 10' = addmargins(table('Overlaps sig cut 10' = countOverlaps(x[[1]], sig[[2]]) > 0, 'Sig cut 5' = x[[1]]$significantFWER == 'TRUE')),
        '10 versus 5' = addmargins(table('Overlaps sig cut 5' = countOverlaps(x[[2]], sig[[1]]) > 0, 'Sig cut 10' = x[[2]]$significantFWER == 'TRUE'))
    )
    percent <- lapply(counts, function(z) { round(z / z[3, 3] * 100, 2) })
    return(list('counts' = counts, 'percent' = percent))
})

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()