#!/bin/sh

## Usage
# sh step5-derfinderReport.sh shulha run8-v1.5.35
# sh step5-derfinderReport.sh epimap run1-v1.5.38 H3K27ac
# sh step5-derfinderReport.sh epimap run1-v1.5.38 H3K4me3

# Define variables
EXPERIMENT=$1
SHORT="derR-${EXPERIMENT}"
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
#$ -l mem_free=200G,h_vmem=250G,h_fsize=20G
#$ -N ${sname}
#$ -hold_jid derM-${EXPERIMENT}.${PREFIX}*

echo "**** Job starts ****"
date

mkdir -p ${WDIR}/${outdir}/logs

# merge results
cd ${WDIR}
module load R/3.3
Rscript -e "library(regionReport); load('${MAINDIR}/CoverageInfo/fullCov.Rdata'); derfinderReport(prefix='${outdir}', browse=FALSE, nBestRegions = 100,  nBestClusters = 20, fullCov=fullCov, device='CairoPNG', clean = FALSE); Sys.time(); proc.time(); options(width = 120); devtools::session_info()"

# Move log files into the logs directory
mv ${ROOTDIR}/${sname}.* ${WDIR}/${outdir}/logs/

echo "**** Job ends ****"
date
EOF
call="qsub .${sname}.sh"
echo $call
$call
