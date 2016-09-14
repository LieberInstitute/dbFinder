## Usage:
# Rscript make_pd.R
# mkdir -p logs
# Rscript make_pd.R > logs/make_pd_log.txt 2>&1

## Libraries needed
library('devtools')

message(paste(Sys.time(), 'loading data'))
## Read log file containing the total mapped reads info
logs <- readLines('logs/makeBai-USC.o5773966')
## Drop job start/end and timing lines
logs <- logs[c(3:(length(logs) - 2))]
mapped <- as.integer(logs[rep(c(FALSE, TRUE), length(logs) / 2)])
names(mapped) <- logs[rep(c(TRUE, FALSE), length(logs) / 2)]
names(mapped) <- gsub('.*/|.bam', '', names(mapped))

## Read clinical data
clinical <- read.csv('/dcl01/lieber/ajaffe/psychENCODE_Data/USC_U01MH103346/USC_CNON-U01MH103346_Clinical_Metadata_August2016Release.csv')

## Read ChIP-seq data
chip <- read.csv('/dcl01/lieber/ajaffe/psychENCODE_Data/USC_U01MH103346/USC_CNON-U01MH103346_ChIP-seq_Metadata_August2016Release-final.csv')


message(paste(Sys.time(), 'merging data'))
pd <- merge(chip, clinical, by = 'Individual_ID')
map <- match(pd$File_Name, names(mapped))

print("Samples whose file names don't exactly match")
pd$File_Name[is.na(map)]
names(mapped)[!seq_len(length(mapped)) %in% map]

## Match with some replacements
map[is.na(map)] <- match(gsub('ChIPseq_', '', gsub('Epigenomics', 'RNA', pd$File_Name[is.na(map)])), names(mapped))

## Sample completely missing
print('Missing sample')
pd[is.na(map), ]

## Add mapped info
pd$TotalMapped <- mapped[map]

## Number of million of mapped reads
print('Number of mapped reads in millions')
summary(pd$TotalMapped / 1e6)

unmapped <- pd$TotalReads - pd$TotalMapped
## Number of million of unmapped reads
print('Number of unmapped reads in millions')
summary(unmapped / 1e6)

## Weird sample: note TotalReads is < UniquelyMappedReads in this case
print('Weirdly annotated sample: probably an error in TotalReads')
pd[unmapped < 0, ]

## Save data
message(paste(Sys.time(), 'saving pd.Rdata object'))
save(pd, file = 'pd.Rdata')

## Reproducibility info
proc.time()
options(width = 120)
session_info()