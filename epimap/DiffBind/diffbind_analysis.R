## Run DiffBind analysis

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

## Select input samples
cell_types <- c('NeuN-', 'NeuN+')
inputs <- lapply(cell_types, function(cell) {
    subset(pd, CellType == cell & HistoneMark == 'Input')
})
names(inputs) <- cell_types

input_bam <- lapply(inputs, function(x) x$bamFile)
input_id <- lapply(inputs, function(x) x$Sample_ID)

pd <- pd[colsubset, ]
stopifnot(all(pd$HistoneMark == opt$histone))

## Change names to avoid issues later on
cond <- ifelse(pd$CellType == 'NeuN-', 'NeunFalse', 'NeunTrue')
fac <- ifelse(substr(gsub('.*_', '', groupInfo), 7, 7) == ')', 'Age1', 'Age2')

## Build sampleSheet
sampleSheet <- data.frame(
    SampleID = pd$Sample_ID,
    Tissue = pd$BrainRegion,
    Factor = fac,
    Condition = cond,
    Treatment = paste0(pd$BrainRegion, cond, fac),
    bamReads = pd$bamFile,
    bamControl = unlist(input_bam[pd$CellType], use.names = FALSE),
    ControlID = unlist(input_id[pd$CellType], use.names = FALSE),
    Peaks = file.path(macsdir, paste0(pd$Sample_ID, '_macs_out_peaks.xls')),
    PeakCaller = 'macs',
    PeakFormat =  'macs', stringsAsFactors = FALSE
)

## For testing
if(FALSE) sampleSheet <- head(sampleSheet)

message(paste(Sys.time(), 'creating DBA object'))
minOv <- 2
diffbind <- dba(sampleSheet = sampleSheet, minOverlap = minOv)
diffbind


## Without using summit
message(paste(Sys.time(), 'counting without specifying summits'))
system.time( diffbind_summit_no <- dba.count(diffbind, minOverlap = minOv) )

message(paste(Sys.time(), 'saving results (no summits) -- in case code breaks'))
save(diffbind_summit_no, file = file.path(maindir, paste0(opt$histone,
    '_diffbind_summit_no.Rdata')))

message(paste(Sys.time(), 'setting contrast'))
diffbind_summit_no <- dba.contrast(diffbind_summit_no, categories=DBA_CONDITION)

message(paste(Sys.time(), 'performing analysis'))
system.time( diffbind_summit_no <- dba.analyze(diffbind_summit_no) )

print('diffbind_summit_no size information')
print(object.size(diffbind_summit_no), units = 'Mb')

message(paste(Sys.time(), 'saving results (no summits)'))
save(diffbind_summit_no, file = file.path(maindir, paste0(opt$histone,
    '_diffbind_summit_no.Rdata')))

message(paste(Sys.time(), 'extracting results'))
diffbind_summit_no_res <- dba.report(diffbind_summit_no, th = 1)
diffbind_summit_no_res$FWER <- p.adjust(mcols(diffbind_summit_no_res)[['p-value']],
    method = 'bonferroni')

message(paste(Sys.time(), 'saving results report (no summits)'))
save(diffbind_summit_no_res, file = file.path(maindir, paste0(opt$histone,
    '_diffbind_summit_no_res.Rdata')))

print('Number of dbPeaks with FWER < 0.05')
table(diffbind_summit_no_res$FWER < 0.05)
round(table(diffbind_summit_no_res$FWER < 0.05) /
    length(diffbind_summit_no_res) * 100, 2)

print('Width summary for all peaks')
summary(width(diffbind_summit_no_res))

print('Width summary for non-sig peaks')
summary(width(diffbind_summit_no_res[diffbind_summit_no_res$FWER >= 0.05]))

print('Width summary for sig peaks')
summary(width(diffbind_summit_no_res[diffbind_summit_no_res$FWER < 0.05]))


## Using summits
message(paste(Sys.time(), 'counting with summits'))
system.time( diffbind_summit <- dba.count(diffbind, minOverlap = minOv, summits = 250) )

message(paste(Sys.time(), 'setting contrast'))
diffbind_summit <- dba.contrast(diffbind_summit, categories=DBA_CONDITION)

message(paste(Sys.time(), 'performing analysis'))
diffbind_summit <- dba.analyze(diffbind_summit)

print('diffbind_summit size information')
print(object.size(diffbind_summit), units = 'Mb')

message(paste(Sys.time(), 'saving results (summits)'))
save(diffbind_summit, file = file.path(maindir, paste0(opt$histone,
    '_diffbind_summit.Rdata')))

message(paste(Sys.time(), 'extracting results'))
diffbind_summit_res <- dba.report(diffbind_summit, th = 1)
diffbind_summit_res$FWER <- p.adjust(mcols(diffbind_summit_res)[['p-value']],
    method = 'bonferroni')
    
message(paste(Sys.time(), 'saving results report (summits)'))
save(diffbind_summit_res, file = file.path(maindir, paste0(opt$histone,
    '_diffbind_summit_res.Rdata')))
    
print('Number of dbPeaks with FWER < 0.05')
table(diffbind_summit_res$FWER < 0.05)
round(table(diffbind_summit_res$FWER < 0.05) / length(diffbind_summit_res) *
    100, 2)

print('Width summary for all peaks')
summary(width(diffbind_summit_res))

print('Width summary for non-sig peaks')
summary(width(diffbind_summit_res[diffbind_summit_res$FWER >= 0.05]))

print('Width summary for sig peaks')
summary(width(diffbind_summit_res[diffbind_summit_res$FWER < 0.05]))

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
