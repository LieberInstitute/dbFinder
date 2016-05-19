## Compare results

## Load libraries
library('getopt')

library('DiffBind')
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

## For testing
if(FALSE) {
    opt <- list(histone = 'H3K4me3')
    opt <- list(histone = 'H3K27ac')
}


## Dirs
rootdir <- '/dcl01/lieber/ajaffe/derRuns/derChIP/epimap/'
derdir <- file.path(rootdir, 'derAnalysis', paste0('run1-v1.5.38-',
    opt$histone))
maindir <- file.path(rootdir, 'DiffBind')
macsdir <- file.path(rootdir, 'macs')
#plotdir <- file.path(maindir, 'plots')
#dir.create(plotdir, showWarnings = FALSE)

## Load sample info
load(file.path(derdir, 'groupInfo.Rdata'))
load(file.path(derdir, 'colsubset.Rdata'))

## Phenotype information
message(paste(Sys.time(), 'loading phenotype information'))
load('/dcl01/lieber/ajaffe/psychENCODE_Data/EpiMap/annotated_phenotype_EpiMap_ChIPseq.rda')
pd <- pd[colsubset, ]
stopifnot(all(pd$HistoneMark == opt$histone))

## Load regions
message(paste(Sys.time(), 'loading fullRegions.Rdata'))
load(file.path(derdir, 'fullRegions.Rdata'))

#keepIndex <- width(fullRegions) >=6
keepIndex <- fullRegions$significantFWER == 'TRUE'
regions <- fullRegions[keepIndex]

## Load DiffBind results
message(paste(Sys.time(), 'loading DiffBind results'))
load(file.path(maindir, paste0(opt$histone, '_diffbind_summit_no_res.Rdata')))
load(file.path(maindir, paste0(opt$histone, '_diffbind_summit_res.Rdata')))

## Group into a list
message(paste(Sys.time(), 'comparing results'))
db <- list('no_summit' = diffbind_summit_no_res, 'summit' = diffbind_summit_res)
print('Number of merged peaks (DiffBind)')
sapply(db, length)

## Subset only significant ones
db_sig <- lapply(db, function(d) { d[d$FWER < 0.5] })
print('Number of significant merged peaks (DiffBind)')
sapply(db_sig, length)

## Overlap between DiffBind results
print('Overlap between DiffBind results (sig only)')
table('no_summit' = countOverlaps(db_sig[[1]], db_sig[[2]]) > 0)
round(table('no_summit' = countOverlaps(db_sig[[1]], db_sig[[2]]) > 0) / length(db_sig[[1]]) * 100, 2)
table('summit' = countOverlaps(db_sig[[2]], db_sig[[1]]) > 0)
round(table('summit' = countOverlaps(db_sig[[2]], db_sig[[1]]) > 0) / length(db_sig[[2]]) * 100, 2)


## Simple overlaps
db_der <- sapply(db_sig, function(d) {
    table(countOverlaps(d, regions) > 0)
})
print('Overlap, query: DiffBind sig')
db_der
round(db_der / colSums(db_der) * 100, 2)

der_db <- sapply(db_sig, function(d) {
    table(countOverlaps(regions, d) > 0)
})
print('Overlap, query: sig dbPeaks')
der_db
round(der_db / colSums(der_db) * 100, 2)


db_der_all <- sapply(db, function(d) {
    table(countOverlaps(d, fullRegions) > 0)
})
print('Overlap, query: DiffBind all, subject: dbPeaks all')
db_der_all
round(db_der_all / colSums(db_der_all) * 100, 2)

der_db_all <- sapply(db, function(d) {
    table(countOverlaps(fullRegions, d) > 0)
})
print('Overlap, query: dbPeaks all, subject: DiffBind all')
der_db_all
round(der_db_all / colSums(der_db_all) * 100, 2)

## Get genome length
data(hg19Ideogram, package = 'biovizBase')
seqlengths(regions) <- seqlengths(hg19Ideogram)[names(seqlengths(regions))]

## Define 1kb tiles
tiles <- tileGenome(seqlengths(regions), tilewidth = 1e3)

sig <- c(list('dbPeaks' = regions), db_sig)
tile_ov <- sapply(sig, function(d) {
    table(countOverlaps(tiles, d) > 0)
})
print('Overlap between 1kb tiles and sig results')
tile_ov
round(tile_ov / colSums(tile_ov) * 100, 2)

all <- c(list('dbPeaks' = fullRegions), db)
tile_ov_all <- sapply(all, function(d) {
    table(countOverlaps(tiles, d) > 0)
})
print('Overlap between 1kb tiles and all results')
tile_ov_all
round(tile_ov_all / colSums(tile_ov_all) * 100, 2)

## Save results
save(db_der, der_db, db_der_all, der_db_all, tile_ov, tile_ov_all,
    file = file.path(maindir, paste0(opt$histone, '_comparison.Rdata')))

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
