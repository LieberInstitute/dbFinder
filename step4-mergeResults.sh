#!/bin/sh

## Usage
# sh step4-mergeResults.sh shulha run8-v1.5.35
# sh step4-mergeResults.sh epimap run1-v1.5.38 H3K27ac
# sh step4-mergeResults.sh epimap run1-v1.5.38 H3K4me3

# Define variables
EXPERIMENT=$1
SHORT="derM-${EXPERIMENT}"
PREFIX=$2
HISTONE=$3

# Directories
ROOTDIR=/dcl01/lieber/ajaffe/derRuns/derChIP
MAINDIR=${ROOTDIR}/${EXPERIMENT}
WDIR=${MAINDIR}/derAnalysis

# Construct shell files
if [[ "${EXPERIMENT}" == "shulha" ]]
then
    outdir="${PREFIX}"
    sname="${SHORT}.${PREFIX}"
elif [[ "${EXPERIMENT}" == 'epimap' ]]
then
    outdir="${PREFIX}-${HISTONE}"
    sname="${SHORT}.${PREFIX}-${HISTONE}"
else
    echo "Specify a valid experiment: shulha, epimap"
fi

echo "Creating script ${sname}"
cat > ${ROOTDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=100G,h_vmem=150G,h_fsize=20G
#$ -N ${sname}
#$ -hold_jid derA-${EXPERIMENT}.${PREFIX}*,fix_ann_chr*

echo "**** Job starts ****"
date

mkdir -p ${WDIR}/${outdir}/logs

# merge results
cd ${WDIR}
module load R/3.3
Rscript -e "library(derfinder); load('/dcl01/lieber/ajaffe/derRuns/derfinderExample/derGenomicState/GenomicState.Hsapiens.UCSC.hg19.knownGene.Rdata'); load('${WDIR}/${outdir}/chr22/optionsStats.Rdata'); chrs <- paste0('chr', c(1:22, 'X', 'Y', 'M')); mergeResults(chrs = chrs, prefix = '${outdir}', genomicState = GenomicState.Hsapiens.UCSC.hg19.knownGene[['fullGenome']], optionsStats = optionsStats); prefix <- '${outdir}'; source('/dcl01/lieber/ajaffe/derRuns/derChIP/fix-colnames.R'); Sys.time(); proc.time(); options(width = 120); devtools::session_info()"

# Move log files into the logs directory
mv ${ROOTDIR}/${sname}.* ${WDIR}/${outdir}/logs/

echo "**** Job ends ****"
date
EOF
call="qsub .${sname}.sh"
echo $call
$call
