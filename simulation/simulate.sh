#!/bin/sh

## Usage
# sh simulate.sh tf
# sh simulate.sh hist


# Define variables
EXPERIMENT=$1
SHORT="sim-${EXPERIMENT}"

# Directories
WDIR=/dcl01/lieber/ajaffe/derRuns/derChIP/simulation

# Construct shell file
echo 'Creating script for loading the Coverage data'
cat > ${WDIR}/.${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=150G,h_vmem=170G,h_fsize=40G
#$ -N ${SHORT}

echo '**** Job starts ****'
date

## Make logs directory
mkdir -p ${WDIR}/logs

## Load required modules
module load macs/2.1.0
module load R/3.3.x
module load samtools/1.1

## Run the simulation
Rscript autosim.R -t "${EXPERIMENT}"

## Move log files into the logs directory
mv ${SHORT}.* logs/

echo '**** Job ends ****'
date
EOF

call="qsub .${SHORT}.sh"
echo $call
$call
