#!/bin/bash	
#$ -cwd
#$ -m e
#$ -l mem_free=2G,h_vmem=4G
#$ -pe local 10
#$ -N fullCov-shula
echo "**** Job starts ****"
date

mkdir -p logs

# Generate HTML
Rscript import-data.R

mv fullCov-shula.* logs/

echo "**** Job ends ****"
date
