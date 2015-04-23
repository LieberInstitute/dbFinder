library('ggplot2')
library('getopt')

## Specify parameters
spec <- matrix(c(
	'experiment', 'e', 1, 'character', 'Experiment',
	'run', 'r', 1, 'character', 'run name',
    'npermute', 'n', 1, 'integer', 'Number of permutations',
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

chrs <- paste0('chr', c(1:22, 'X', 'Y', 'M'))
study <- opt$experiment
run <- opt$run

timediff <- lapply(chrs, function(chr) {
    info <- tryCatch(system(paste0('grep permutation *', study, '*', run, '*', chr, '.e*'), intern = TRUE), warning = function(w) { 'no data'})
    if(info[1] == 'no data') {
        info <- tryCatch(system(paste0('grep permutation ', file.path(study, 'derAnalysis', run, chr, 'logs'), '/*', chr, '.e*'), intern = TRUE), warning = function(w) { 'no data'})
    }
    if(info[1] == 'no data') return(NULL)
    
    time <- strptime(gsub('([[:space:]]*calculate.*$)', '', info),
        format = '%Y-%m-%d %H:%M:%S')
    
    idx <- seq_len(length(info) - 1)
    difftime(time[idx + 1], time[idx], units = 'mins')
})
names(timediff) <- chrs


## Organize time information
chrnum <- gsub('chr', '', chrs)
df <- data.frame(chr = factor(chrnum, levels = chrnum), mean = sapply(timediff, mean), sd = sapply(timediff, sd))

## Group by number of rounds per permutation given the number of chunks & cores used
if(!file.exists(file.path(study, 'derAnalysis', run, 'nChunks.Rdata'))) {
    nChunks <- sapply(chrs, function(chr) { 
        if(!file.exists(file.path(study, 'derAnalysis', run, chr, 'coveragePrep.Rdata')))
            return(NA)
        load(file.path(study, 'derAnalysis', run, chr, 'coveragePrep.Rdata'))
        max(prep$mclapply) 
    })
    save(nChunks, file = file.path(study, 'derAnalysis', run, 'nChunks.Rdata'))
} else {
    load(file.path(study, 'derAnalysis', run, 'nChunks.Rdata'))
}

if(study == 'shula') {
    nCores <- rep(2, 25)
}
names(nCores) <- chrs

df$n <- sapply(timediff, length)
df$se <- df$sd / sqrt(df$n)
df$nChunks <- nChunks
df$nCores <- nCores
df$nRound <- factor(ceiling(nChunks / nCores))


## Print info
rownames(df) <- NULL
print(df)


## Make plot
pdf(file.path(study, 'derAnalysis', run, paste0('permuteTime-', study, '-', run, '.pdf')))
ggplot(df, aes(x = chr, y = mean, color = nRound)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1) + geom_line() + geom_point() + ylab('Time per permutation (minutes)\nMean +- SE') + xlab('Chromosome') + ggtitle(paste('Time info for', study, run)) # + scale_y_continuous(breaks=seq(0, ceiling(max(df$mean + df$se, na.rm = TRUE)), 1))
dev.off()


print('Expected total number of hours per chr and hours remaining')
hours <- data.frame(chr = chrnum, total = round(df$mean * (opt$npermute + 1) / 60, 1), remaining = round(df$mean * (opt$npermute + 1 - df$n - 2 ) / 60, 1))
rownames(hours) <- NULL
print(hours)

print('Expected total number of days per chr and days remaining')
days <- data.frame(chr = chrnum, total = round(df$mean * (opt$npermute + 1) / 60 / 24, 1), remaining = round(df$mean * (opt$npermute + 1 - df$n - 2 ) / 60 / 24, 1))
rownames(days) <- NULL
print(days)




