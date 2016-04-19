#!/bin/sh

## Usage
# sh step3-analyzeChr.sh shulha run6-v1.5.34

# Define variables
EXPERIMENT=$1
SHORT="derA-${EXPERIMENT}"
PREFIX=$2

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
    
    if [[ ${EXPERIMENT} == "shulha" ]]
    then
        CORES=2
    else
        echo "Specify a valid experiment: shulha"
    fi
    
	chr="chr${chrnum}"
	outdir="${PREFIX}/${chr}"
	sname="${SHORT}.${PREFIX}.${chr}"
	cat > ${ROOTDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=8G,h_vmem=10G,h_fsize=10G,h=!compute-04[3-5]*
#$ -N ${sname}
#$ -pe local ${CORES}
#$ -hold_jid derMod-${EXPERIMENT}.${PREFIX}

echo "**** Job starts ****"
date

# Create output directory 
mkdir -p ${WDIR}/${outdir}
# Make logs directory
mkdir -p ${WDIR}/${outdir}/logs

# run analyzeChr()
cd ${WDIR}/${PREFIX}/
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
