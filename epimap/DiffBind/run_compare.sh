#!/bin/sh

## Usage
# sh run_compare.sh

# Directories
ROOTDIR=/dcl01/lieber/ajaffe/derRuns/derChIP/epimap/DiffBind


# Construct shell files
for HISTONE in H3K27ac H3K4me3
    do
    sname="${HISTONE}.compare"
    echo "Creating script ${sname}"
    cat > ${ROOTDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=100G,h_vmem=120G,h_fsize=10G
#$ -N ${sname}

echo "**** Job starts ****"
date

mkdir -p ${ROOTDIR}/logs/

# Compare results
module load R/3.3
Rscript ${ROOTDIR}/compare_results.R -i "${HISTONE}"

# Move log files into the logs directory
mv ${ROOTDIR}/${sname}.* ${ROOTDIR}/logs/

echo "**** Job ends ****"
date
EOF
    ## Run job
    call="qsub .${sname}.sh"
    echo $call
    $call
done
