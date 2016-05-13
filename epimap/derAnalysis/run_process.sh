#!/bin/sh

## Usage
# sh run_process.sh

# Directories
ROOTDIR=/dcl01/lieber/ajaffe/derRuns/derChIP/epimap/derAnalysis


# Construct shell files
for HISTONE in H3K27ac H3K4me3
    do
    if [[ "${HISTONE}" == "H3K27ac" ]]
    then
        CHRNUMS="Y X 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
        
    elif [[ "${HISTONE}" == "H3K4me3" ]]
    then
        CHRNUMS="all"
    fi
    
    for chrnum in ${CHRNUMS}
    do
        chr="chr${chrnum}"
        sname="${HISTONE}.${chr}.process"
        outdir="run1-v1.5.38-${HISTONE}"
        echo "Creating script ${sname}"
        cat > ${ROOTDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=150G,h_vmem=200G,h_fsize=100G
#$ -N ${sname}

echo "**** Job starts ****"
date

mkdir -p ${ROOTDIR}/${outdir}/logs/

# Process results
module load R/3.3
Rscript ${ROOTDIR}/process_results.R -i "${HISTONE}" -c "${chr}"

# Move log files into the logs directory
mv ${ROOTDIR}/${sname}.* ${ROOTDIR}/${outdir}/logs/

echo "**** Job ends ****"
date
EOF
        ## Run job
        call="qsub .${sname}.sh"
        echo $call
        $call
    done
    
    if [[ "${HISTONE}" == "H3K27ac" ]]
    then
        sname="${HISTONE}.merge.process"
        echo "Creating script ${sname}"
        cat > ${ROOTDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=200G,h_vmem=270G,h_fsize=100G
#$ -N ${sname}

echo "**** Job starts ****"
date

mkdir -p ${ROOTDIR}/${outdir}/logs/

# Merge process results
module load R/3.3
Rscript ${ROOTDIR}/merge_process.R

# Move log files into the logs directory
mv ${ROOTDIR}/${sname}.* ${ROOTDIR}/${outdir}/logs/

echo "**** Job ends ****"
date
EOF
        ## Run job
        call="qsub .${sname}.sh"
        echo $call
        $call
    fi
done
