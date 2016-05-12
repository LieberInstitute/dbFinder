#!/bin/sh

## Usage
# sh run_process.sh

# Directories
ROOTDIR=/dcl01/lieber/ajaffe/derRuns/derChIP/epimap/derAnalysis


# Construct shell files
for HISTONE in H3K27ac H3K4me3
    do
    sname="${HISTONE}.process"
    outdir="run1-v1.5.38-${HISTONE}"
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

# Process results
module load R/3.3
Rscript ${ROOTDIR}/process_results.R -i "${HISTONE}"

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
