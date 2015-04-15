#!/bin/bash	
#$ -cwd
#$ -m e
#$ -l mem_free=10G,h_vmem=20G
#$ -pe local 10
#$ -N shula-fullCov
echo "**** Job starts ****"
date

mkdir -p logs

# Generate HTML
Rscript import-data.R

mv shula-fullCov.* logs/

echo "**** Job ends ****"
date
