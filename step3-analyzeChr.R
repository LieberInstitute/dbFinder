## Run derfinder's analysis steps with timing info

## Load libraries
library('getopt')

## Available at http://www.bioconductor.org/packages/release/bioc/html/derfinder.html
library('derfinder')
library('devtools')
library('GenomeInfoDb')

## Specify parameters
spec <- matrix(c(
    'experiment', 'e', 1, 'character', 'Experiment. Either stem, brainspan, snyder, hippo, or simulation',
	'CovFile', 'd', 1, 'character', 'path to the .Rdata file with the results from loadCoverage() with cutoff >= 0',
	'chr', 'c', 1, 'character', 'Chromosome under analysis',
	'mcores', 'm', 1, 'integer', 'Number of cores',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)


## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## Check experiment input
stopifnot(opt$experiment %in% c('shula'))

## Format chromosome name appropriately
if(opt$chr == 'M') opt$chr <- 'chrM' else opt$chr <- mapSeqlevels(opt$chr, 'UCSC')


message('Loading Rdata file with the output from loadCoverage()')
load(opt$CovFile)

## Make it easy to use the name later. Here I'm assuming the names were generated using output='auto' in loadCoverage()
eval(parse(text=paste0('covData <- ', opt$chr, 'CovInfo')))
eval(parse(text=paste0('rm(', opt$chr, 'CovInfo)')))


## Load the models
load('models.Rdata')

## Load group information
load('groupInfo.Rdata')

if(file.exists('colsubset.Rdata')) {
    load('colsubset.Rdata')
} else {
    colsubset <- NULL
}


## Run the analysis
if(opt$experiment == 'shula') {
    analyzeChr(chr = opt$chr, coverageInfo = covData, models = models, 
        cutoffFstat = 0.01, colsubset = colsubset, cutoffPre = 2,
        nPermute = 100, seeds = seq_len(100) + 150422, maxClusterGap = 3000,
        groupInfo = groupInfo, mc.cores = opt$mcores,
        lowMemDir = file.path(tempdir(), opt$chr, 'chunksDir'))    
}


## Done
proc.time()
options(width = 120)
session_info()
