#!/bin/bash	
#$ -cwd
#$ -m e
#$ -l mem_free=10G,h_vmem=20G
#$ -pe local 10
#$ -N shula-fullCov
echo "**** Job starts ****"
date

# Generate HTML
Rscript import-data.R

echo "**** Job ends ****"
date
