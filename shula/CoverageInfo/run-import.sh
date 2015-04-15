#!/bin/bash	
#$ -cwd
#$ -m e
#$ -l mem_free=10G,h_vmem=20G
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
