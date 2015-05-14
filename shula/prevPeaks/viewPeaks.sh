#!/bin/bash	
#$ -cwd
#$ -m e
#$ -l mem_free=30G,h_vmem=100G
#$ -N viewPeaks
echo "**** Job starts ****"
date

# Generate HTML
module load R/3.1.x
Rscript -e "library(rmarkdown); render('viewPeaks', clean = FALSE)"

echo "**** Job ends ****"
date
