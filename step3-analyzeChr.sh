#!/bin/sh

## Usage
# sh step3-analyzeChr.sh shulha run8-v1.5.35
# sh step3-analyzeChr.sh epimap run1-v1.5.38 H3K27ac
# sh step3-analyzeChr.sh epimap run1-v1.5.38 H3K4me3

# Define variables
EXPERIMENT=$1
SHORT="derA-${EXPERIMENT}"
PREFIX=$2
HISTONE=$3

# Directories
ROOTDIR=/dcl01/lieber/ajaffe/derRuns/derChIP
MAINDIR=${ROOTDIR}/${EXPERIMENT}
WDIR=${MAINDIR}/derAnalysis
DATADIR=${MAINDIR}/CoverageInfo

# Construct shell files
CHRNUMS="M 22 21 Y 20 19 18 17 16 15 14 13 12 11 10 9 8 X 7 6 5 4 3 2 1"

for chrnum in ${CHRNUMS}
do
	echo "Creating script for chromosome ${chrnum}"
    chr="chr${chrnum}"
    
    if [[ "${EXPERIMENT}" == "shulha" ]]
    then
        CORES=4
    	outdir="${PREFIX}/${chr}"
        rundir="${PREFIX}"
    	sname="${SHORT}.${PREFIX}.${chr}"
    elif [[ "${EXPERIMENT}" == "epimap" ]]
    then
        CORES=20
    	outdir="${PREFIX}-${HISTONE}/${chr}"
        rundir="${PREFIX}-${HISTONE}"
    	sname="${SHORT}.${PREFIX}-${HISTONE}.${chr}"
    else
        echo "Specify a valid experiment: shulha, epimap"
    fi
	
	cat > ${ROOTDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=2G,h_vmem=3G,h_fsize=10G,h=!compute-04[3-5]*
#$ -N ${sname}
#$ -pe local ${CORES}
#$ -hold_jid derMod-${EXPERIMENT}.${PREFIX}*

echo "**** Job starts ****"
date

# Create output directory 
mkdir -p ${WDIR}/${outdir}
# Make logs directory
mkdir -p ${WDIR}/${outdir}/logs

# run analyzeChr()
cd ${WDIR}/${rundir}/
module load R/3.3
Rscript ${ROOTDIR}/step3-analyzeChr.R -d "${DATADIR}/${chr}CovInfo-filtered.Rdata" -c "${chrnum}" -m ${CORES} -e "${EXPERIMENT}"

# Move log files into the logs directory
mv ${ROOTDIR}/${sname}.* ${WDIR}/${outdir}/logs/

echo "**** Job ends ****"
date
EOF
	call="qsub .${sname}.sh"
	$call
done
