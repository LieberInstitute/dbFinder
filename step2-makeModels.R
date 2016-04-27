## Load the data without a filter, save it, then filter it for derfinder processing steps

## Load libraries
library('getopt')

## Available at http://www.bioconductor.org/packages/release/bioc/html/derfinder.html
library('derfinder')
library('devtools')

## Specify parameters
spec <- matrix(c(
	'experiment', 'e', 1, 'character', 'Experiment. Either shulha or epimap',
    'histone', 'i', 2, 'character', 'For epimap, the histone mark to use. Either H3K27ac or H3K4me3.',
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
stopifnot(opt$experiment %in% c('shulha', 'epimap'))

if(opt$experiment != 'epimap') {
    ## Load the coverage information
    load(file.path('..', '..', 'CoverageInfo', 'fullCov.Rdata'))
    load(file.path('..', '..', 'CoverageInfo', 'chr22CovInfo.Rdata'))

    ## Identify the samplefiles
    files <- colnames(chr22CovInfo$coverage)
}

 ## Calculate the library adjustments and build the models
buildModels <- function(fullCov, testvars, colsubset = NULL, adjustvars = NULL) {
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
        adjustvars = adjustvars, testIntercept = FALSE)
    
    return(models)
}


if(opt$experiment == 'shulha') {

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
} else if (opt$experiment == 'epimap') {
   ## Load the phenotype data
    load('/dcl01/lieber/ajaffe/psychENCODE_Data/EpiMap/annotated_phenotype_EpiMap_ChIPseq.rda')
    
    ## For testing
    if(FALSE) opt <- list(histone = 'H3K27ac')
    
    stopifnot(opt$histone %in% c('H3K4me3', 'H3K27ac'))
    ## Subset to apropriate mark
    colsubset <- which(pd$HistoneMark == opt$histone)
    
    ## Save colsubset
    save(colsubset, file = 'colsubset.Rdata')
    
    ## Define test variables
    testvars <- pd[colsubset, c('BrainRegion', 'CellType', 'AgeDeath')]
    
    ## Define models
    models <- list(
        mod = model.matrix(~ testvars$BrainRegion + testvars$AgeDeath + testvars$CellType),
        mod0 = model.matrix(~0 + rep(1, nrow(testvars)))
    )
    
    ## Define groups
    groupvars <- testvars
    groupvars$AgeDeath <- Hmisc::cut2(testvars$AgeDeath, g = 2)
    groupInfo <- factor(apply(groupvars, 1, function(x) paste(x, collapse = '_')))
}

## Save models
save(models, file="models.Rdata")

## Save information used for analyzeChr(groupInfo)
save(groupInfo, file="groupInfo.Rdata")

## Done :-)
proc.time()
options(width = 120)
session_info()
