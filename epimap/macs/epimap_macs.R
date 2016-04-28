## This script generates shell files that run MACS on the samples
## Usage:
# module load R/3.3
# mkdir -p logs
# Rscript epimap_macs.R > logs/epimap_macs_log.txt 2>&1

## Define path
wdir <- '/dcl01/lieber/ajaffe/derRuns/derChIP/epimap/macs'

## Load phenotype table
load('/dcl01/lieber/ajaffe/psychENCODE_Data/EpiMap/annotated_phenotype_EpiMap_ChIPseq.rda')

## Marks to work with
cell_types <- c('NeuN-', 'NeuN+')

for(cell in cell_types) {
    pd_subset <-  subset(pd, CellType == cell & HistoneMark != 'Input')
    input <- subset(pd, CellType == cell & HistoneMark == 'Input')
    stopifnot(nrow(input) == 1)
    for(sample in pd_subset$Sample_ID) {
        job_name <- paste0('macs-', sample)
        message(paste(Sys.time(), 'creating script for job', job_name))
        script_file <- paste0('.', job_name, '.sh')
        cat(paste0("#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=15G,h_vmem=20G,h_fsize=20G
#$ -N ",job_name, "

echo \"**** Job starts ****\"
date

## Create logs dir
mkdir -p ", wdir, "/logs

# run MACS
cd ", wdir, "
module load macs/2.1.0
macs2 callpeak -t ", pd_subset$bamFile[pd_subset$Sample_ID == sample], " -c ", input$bamFile, " --tsize 100 --bw 230 -n ", sample, "_macs_out

# Move log files into the logs directory
mv ", wdir, "/", job_name, ".* ", wdir, "/logs/

echo \"**** Job ends ****\"
date
"), file = script_file)
        message(paste(Sys.time(), 'qsub', script_file))
        system(paste('qsub', script_file))
    }
}

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
