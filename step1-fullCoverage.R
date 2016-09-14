## Load the data without a filter, save it, then filter it for derfinder processing steps

## Load libraries
library('getopt')

## Available at http://www.bioconductor.org/packages/release/bioc/html/derfinder.html
library('derfinder')
library('BiocParallel')
library('devtools')

## Specify parameters
spec <- matrix(c(
	'datadir', 'd', 1, 'character', 'Data directory, matched with rawFiles(datadir)',
	'pattern', 'p', 1, 'character', 'Sample pattern',
	'cutoff', 'c', 1, 'numeric', 'Filtering cutoff used',
	'mcores', 'm', 1, 'integer', 'Number of cores',
    'fileStyle', 'f', 2, 'character', 'FileStyle used for naming the chromosomes',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## Default values
if (is.null(opt$fileStyle)) opt$fileStyle <- 'UCSC'
if (is.null(opt$totalMapped)) opt$totalMapped <- NULL
if (is.null(opt$targetSize)) opt$targetSize <- 80e6

## Identify the data directories
if(opt$datadir == '/dcs01/ajaffe/ChIPseq/Shulha2013/BED') {
    #files <- rawFiles(datadir=opt$datadir, samplepatt=opt$pattern)
    beds <- dir(opt$datadir, pattern = opt$pattern, full.names = TRUE)
    names(beds) <- dir(opt$datadir, pattern = opt$pattern)
    
    stop('Use shulha/CoverageInfo/run-import.sh instead')
} else if (opt$datadir == '/dcl01/lieber/ajaffe/psychENCODE_Data/EpiMap') {
    load('/dcl01/lieber/ajaffe/psychENCODE_Data/EpiMap/annotated_phenotype_EpiMap_ChIPseq.rda')
    files <- pd$bamFile
    names(files) <- pd$Sample_ID
    opt$totalMapped <- pd$totalMapped
    print(summary(pd$totalMapped) / 1e6)
} else if (opt$datadir == '/dcl01/lieber/ajaffe/psychENCODE_Data/USC_U01MH103346') {
    load('/dcl01/lieber/ajaffe/derRuns/derChIP/USC/pd.Rdata')
    pd <- pd[!is.na(pd$bamFile), ]
    files <- pd$bamFile
    names(files) <- pd$Sample_ID
    opt$totalMapped <- pd$TotalMapped
    print(summary(pd$TotalMapped) / 1e6)
} else if (opt$datadir == '') {
}


## Load the coverage information without filtering
chrs <- paste0('chr', c(1:22, 'X', 'Y', 'M'))

fullCov <- fullCoverage(files = files, chrs = chrs, mc.cores = opt$mcores, fileStyle = opt$fileStyle, outputs = 'auto')

message(paste(Sys.time(), 'Saving the full (unfiltered) coverage data'))
save(fullCov, file='fullCov.Rdata')

rm(fullCov)

fullCov_files <- as.list(dir(pattern = 'chr'))
names(fullCov_files) <- sapply(fullCov_files, function(x) gsub('CovInfo.Rdata', '', x))

## Filter the data and save it by chr
myFilt <- function(chr, rawData_file, cutoff, totalMapped = NULL, targetSize = 80e6) {
    library('derfinder')
    
    ## Load raw data
    message(paste(Sys.time(), 'Loading raw file', rawData_file, 'for', chr))
    load(rawData_file)
    rawData <- get(paste0(chr, 'CovInfo'))
    
	## Filter the data
    message(paste(Sys.time(), 'Filtering chromosome', chr))
	res <- filterData(data = rawData$coverage, cutoff = cutoff, index = NULL,
        totalMapped = totalMapped, targetSize = targetSize)
	
	## Save it in a unified name format
	varname <- paste0(chr, 'CovInfo')
	assign(varname, res)
	output <- paste0(varname, '-filtered.Rdata')
	
	## Save the filtered data
	save(list = varname, file = output, compress='gzip')
	
	## Finish
	return(invisible(NULL))
}


message(paste(Sys.time(), 'Filtering and saving the data with cutoff', opt$cutoff))
filteredCov <- bpmapply(myFilt, names(fullCov_files), fullCov_files, BPPARAM = SnowParam(workers = opt$mcores), MoreArgs = list(cutoff = opt$cutoff, totalMapped = opt$totalMapped, targetSize = opt$targetSize))

## Done!
proc.time()
options(width = 120)
session_info()
