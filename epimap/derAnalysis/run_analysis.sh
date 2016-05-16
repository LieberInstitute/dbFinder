#!/bin/sh

## Usage
# sh run_analysis.sh

# Directories
ROOTDIR=/dcl01/lieber/ajaffe/derRuns/derChIP/epimap/derAnalysis


# Construct shell files
for HISTONE in H3K27ac H3K4me3
    do
    sname="${HISTONE}.analysis"
    outdir="run1-v1.5.38-${HISTONE}"
    echo "Creating script ${sname}"
    cat > ${ROOTDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=20G,h_vmem=30G,h_fsize=100G
#$ -pe local 10
#$ -N ${sname}
#$ -hold_jid ${HISTONE}.merge.process

echo "**** Job starts ****"
date

mkdir -p ${ROOTDIR}/${outdir}/logs/

# Analyze results
module load R/3.3
Rscript ${ROOTDIR}/analyze_results.R -i "${HISTONE}"

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
