#!/bin/bash	
#$ -cwd
#$ -m e
#$ -l mem_free=30G,h_vmem=100G
#$ -N visPeaks
#$ -hold_jid derM-shulha
echo "**** Job starts ****"
date

mkdir -p logs

# Generate HTML
module load R/3.3
Rscript -e "library(rmarkdown); render('viewPeaks.Rmd', clean = FALSE)"

mv visPeaks.* logs/

echo "**** Job ends ****"
date
