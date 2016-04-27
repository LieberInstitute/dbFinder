#!/bin/sh

## Usage
# sh step2-makeModels.sh shulha run8-v1.5.35
# sh step2-makeModels.sh epimap run1-v1.5.38 H3K27ac
# sh step2-makeModels.sh epimap run1-v1.5.38 H3K4me3

# Define variables
EXPERIMENT=$1
SHORT="derMod-${EXPERIMENT}"
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
#$ -l mem_free=50G,h_vmem=100G,h_fsize=10G
#$ -N ${sname}
#$ -hold_jid fullCov-${EXPERIMENT}

echo "**** Job starts ****"
date

mkdir -p ${WDIR}/${outdir}/logs

# merge results
cd ${WDIR}/${outdir}/
module load R/3.3

Rscript ${ROOTDIR}/step2-makeModels.R -e "${EXPERIMENT}" -i "${HISTONE}"

# Move log files into the logs directory
mv ${ROOTDIR}/${sname}.* ${WDIR}/${outdir}/logs/

echo "**** Job ends ****"
date
EOF


## Run job
call="qsub .${sname}.sh"
echo $call
$call
