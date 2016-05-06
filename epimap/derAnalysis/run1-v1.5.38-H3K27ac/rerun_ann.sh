#!/bin/sh

## Usage
# sh rerun_ann.sh chr2
# sh rerun_ann.sh chr5
# sh rerun_ann.sh chr11
# sh rerun_ann.sh chr16

# Define variables
CHR=$1

# Directories
WDIR=/dcl01/lieber/ajaffe/derRuns/derChIP/epimap/derAnalysis/run1-v1.5.38-H3K27ac

## Construct shell files
sname="fix_ann_${CHR}"
echo "Creating script ${sname}"
cat > ${WDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=80G,h_vmem=100G,h_fsize=10G
#$ -N ${sname}

echo "**** Job starts ****"
date

mkdir -p ${WDIR}/${CHR}/logs

## Fix annotation
module load R/3.3
Rscript ${WDIR}/rerun_annotation.R -c "${CHR}"

# Move log files into the logs directory
mv ${WDIR}/${sname}.* ${WDIR}/${CHR}/logs/

echo "**** Job ends ****"
date
EOF


## Run job
call="qsub .${sname}.sh"
echo $call
$call
