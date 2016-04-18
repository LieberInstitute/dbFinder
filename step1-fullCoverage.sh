#!/bin/sh


## Usage
# sh step1-fullCoverage.sh shulha


# Define variables
EXPERIMENT=$1
SHORT="fullCov-${EXPERIMENT}"
CORES=10

# Directories
ROOTDIR=/dcl01/lieber/ajaffe/derRuns/derChIP
MAINDIR=${ROOTDIR}/${EXPERIMENT}
WDIR=${MAINDIR}/CoverageInfo

if [[ "${EXPERIMENT}" == "shulha" ]]
then
    DATADIR=/dcs01/ajaffe/ChIPseq/Shulha2013/BED
    PATTERN='c'
    CUTOFF=2
else
    echo "Specify a valid experiment: shulha"
fi



# Construct shell file
echo 'Creating script for loading the Coverage data'
cat > ${ROOTDIR}/.${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=10G,h_vmem=15G,h_fsize=40G
#$ -N ${SHORT}
#$ -pe local ${CORES}

echo '**** Job starts ****'
date

# Make logs directory
mkdir -p ${WDIR}/logs

# Load the data, save the coverage without filtering, then save each file separately
cd ${WDIR}
module load R/3.3
Rscript ${ROOTDIR}/step1-fullCoverage.R -d "${DATADIR}" -p "${PATTERN}" -c "${CUTOFF}" -m ${CORES}

## Move log files into the logs directory
mv ${ROOTDIR}/${SHORT}.* ${WDIR}/logs/

echo '**** Job ends ****'
date
EOF

call="qsub .${SHORT}.sh"
echo $call
$call
