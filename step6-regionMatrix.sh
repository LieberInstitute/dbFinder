#!/bin/sh

## Usage
# sh step6-regionMatrix.sh epimap

# Define variables
EXPERIMENT=$1
SHORT="regMat-${EXPERIMENT}"

# Directories
ROOTDIR=/dcl01/lieber/ajaffe/derRuns/derChIP
MAINDIR=${ROOTDIR}/${EXPERIMENT}
WDIR=${MAINDIR}/regionMatrix

if [[ "${EXPERIMENT}" == "epimap" ]]
then
    CUTOFFS="5 10"
    RLENGTH=100
else
    echo "Specify a valid experiment: epimap"
fi


# Construct shell files


for CUTOFF in ${CUTOFFS}
do
    for HISTONE in H3K27ac H3K4me3
    do
        for chrnum in 22 21 Y 20 19 18 17 16 15 14 13 12 11 10 9 8 X 7 6 5 4 3 2 1
        do
            chr="chr${chrnum}"
            sname="${SHORT}-${HISTONE}-cut-${CUTOFF}-${chr}"
            echo "Creating script ${sname}"

            cat > ${ROOTDIR}/.${sname}.sh <<EOF
#!/bin/bash	
#$ -cwd
#$ -m e
#$ -l mem_free=70G,h_vmem=90G,h_fsize=30G
#$ -N ${sname}
#$ -hold_jid fullCov-${EXPERIMENT}

echo "**** Job starts ****"
date

# Make logs directory
mkdir -p ${WDIR}/logs

# Load coverage & get region matrix
cd ${WDIR}
module load R/3.3
Rscript ${ROOTDIR}/step6-regionMatrix.R -m "${MAINDIR}" -c "${chr}" -r ${RLENGTH} -t ${CUTOFF} -i "${HISTONE}"

## Move log files into the logs directory
mv ${ROOTDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date
EOF

            call="qsub .${sname}.sh"
            echo $call
            $call
        done
    
        sname="${SHORT}-${HISTONE}-cut-${CUTOFF}-merge"
    
        echo "Creating script for merging the regionMatrix results for cutoff ${CUTOFF}"
        cat > ${ROOTDIR}/.${sname}.sh <<EOF
#!/bin/bash	
#$ -cwd
#$ -m e
#$ -l mem_free=100G,h_vmem=120G,h_fsize=40G
#$ -N ${sname}
#$ -hold_jid ${SHORT}-${HISTONE}-cut-${CUTOFF}-chr*

echo "**** Job starts ****"
date

# Make logs directory
mkdir -p ${WDIR}/logs

# Load coverage & get region matrix
cd ${WDIR}
module load R/3.3
Rscript ${ROOTDIR}/step6b-regionMatrix-merge.R -t ${CUTOFF} -i "${HISTONE}"

## Move log files into the logs directory
mv ${ROOTDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date
EOF
    
        call="qsub .${sname}.sh"
        echo $call
        $call
    done
done
