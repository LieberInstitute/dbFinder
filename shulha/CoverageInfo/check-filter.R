setwd('/dcs01/ajaffe/Brain/derRuns/derChIP/shulha/CoverageInfo')

library('GenomicRanges')
chrs <- paste0('chr', c(1:22, 'X', 'Y', 'M'))

## For testing:
# chrs <- c('chrM', 'chrY')

filtered <- lapply(chrs, function(chr) {
    load(paste0(chr, 'CovInfo.Rdata'))
    eval(parse(text=paste0('covData <- ', chr, 'CovInfo')))
    #eval(parse(text=paste0('rm(', chr, 'CovInfo)')))
    data.frame('Remaining' = sum(covData$position), 'Total' = length(covData$position), 'chr' = chr)
})
filtered <- do.call(rbind, filtered)

filtered$RemainPercent <- filtered$Remaining / filtered$Total * 100
filtered$RemovedPercent <- 100 - filtered$RemainPercent
filtered

## Total percent remaining
sum(as.numeric(filtered$Remaining)) / sum(as.numeric(filtered$Total)) * 100
filtered$Remaining / 1e6

save(filtered, file = 'filtered.Rdata')
