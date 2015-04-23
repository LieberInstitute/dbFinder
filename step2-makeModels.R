## Load the data without a filter, save it, then filter it for derfinder processing steps

## Load libraries
library('getopt')

## Available at http://www.bioconductor.org/packages/release/bioc/html/derfinder.html
library('derfinder')
library('devtools')

## Specify parameters
spec <- matrix(c(
	'experiment', 'e', 1, 'character', 'Experiment. Either stem, brainspan, snyder, hippo, or simulation',
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

if(opt$experiment != 'brainspan') {
    ## Load the coverage information
    load(file.path('..', '..', 'CoverageInfo', 'fullCov.Rdata'))
    load(file.path('..', '..', 'CoverageInfo', 'chr22CovInfo.Rdata'))

    ## Identify the samplefiles
    files <- colnames(chr22CovInfo$coverage)
}

 ## Calculate the library adjustments and build the models
buildModels <- function(fullCov, testvars, colsubset = NULL) {
    ## Determine sample size adjustments
    if(file.exists("sampleDepths.Rdata")) {
    	load("sampleDepths.Rdata")
    } else {
    	if(file.exists("collapsedFull.Rdata")) {
    		load("collapsedFull.Rdata")
    	} else {
    		## Collapse
    		collapsedFull <- collapseFullCoverage(fullCov, colsubset = colsubset, save=TRUE)
    	}

    	## Get the adjustments
    	sampleDepths <- sampleDepth(collapsedFull = collapsedFull, probs = 1,
            nonzero = TRUE, scalefac = 1, center = FALSE)
    	save(sampleDepths, file="sampleDepths.Rdata")
    }
    ## Build the models
    models <- makeModels(sampleDepths = sampleDepths, testvars = testvars,
        adjustvars = NULL, testIntercept = FALSE)
    
    return(models)
}


if(opt$experiment == 'shula') {

    ## Load the information table
    pd <- read.csv('/home/epi/ajaffe/Lieber/Projects/ChIP-Seq/chip_phenotype.csv')
    pd$Sample <- gsub('http://zlab.umassmed.edu/zlab/publications/ShulhaPLOSGen2013/|p.zip', '', pd$Filename)
    ## This also works:
    # gsub('put', '', gsub('-', 'N', tolower(pd$Sample.ID)))

    ## Note one name doesn't match, since it has a p left
    pd$Sample[!pd$Sample %in% files]
    files[!files %in% pd$Sample]
    files <- gsub('p', '', files)
    
    ## Reorder pd
    pd <- pd[match(files, pd$Sample), ]
    
    ## Get age
    pd$AgeYr <- as.numeric(gsub(' gw| yr', '', pd$Age))
    pd$AgeYr[grepl('gw', pd$Age)] <- (as.integer(gsub(' gw', '', pd$Age[grepl('gw', pd$Age)])) - 42) / 52
    
    ## Drop input and NeuN- samples
    colsubset <- which(!grepl('N|in', pd$Sample))
    save(colsubset, file = 'colsubset.Rdata')
    
    testvars <- pd$AgeYr[colsubset]
    ## For age with bs-splines
    # library('splines')
    # testvars <- bs(pd$AgeYr[colsubset], df = 5)
    
    ## Define the groups
    groupInfo <- cut(testvars, breaks = c(-1, 0, 1, 10, 20, 30, 100))
    
    ## Build models
    models <- buildModels(fullCov, testvars, colsubset)
}

## Save models
save(models, file="models.Rdata")

## Save information used for analyzeChr(groupInfo)
save(groupInfo, file="groupInfo.Rdata")

## Done :-)
proc.time()
options(width = 120)
session_info()
