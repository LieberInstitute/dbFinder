#!/bin/sh

## Usage
# sh run_analysis.sh

# Directories
ROOTDIR=/dcl01/lieber/ajaffe/derRuns/derChIP/epimap/regionMatrix


# Construct shell files
for HISTONE in H3K27ac H3K4me3
    do
        for CUTOFF in 5 10
        do
        sname="${HISTONE}.cut${CUTOFF}.ER-analysis"
        outdir="${HISTONE}-cut${CUTOFF}"
        echo "Creating script ${sname}"
        cat > ${ROOTDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=4G,h_vmem=5G,h_fsize=100G
#$ -pe local 10
#$ -N ${sname}
#$ -hold_jid regMat-epimap-${HISTONE}-cut-${CUTOFF}-merge

echo "**** Job starts ****"
date

mkdir -p ${ROOTDIR}/${outdir}
mkdir -p ${ROOTDIR}/${outdir}/logs/

# Analyze results
module load R/3.3
Rscript ${ROOTDIR}/analyze_results.R -i "${HISTONE}" -t ${CUTOFF}

# Move log files into the logs directory
mv ${ROOTDIR}/${sname}.* ${ROOTDIR}/${outdir}/logs/

echo "**** Job ends ****"
date
EOF
        ## Run job
        call="qsub .${sname}.sh"
        echo $call
        $call
        
        sname="${HISTONE}.cut${CUTOFF}.ER-analysis.all"
        echo "Creating script ${sname}"
        cat > ${ROOTDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=4G,h_vmem=5G,h_fsize=100G
#$ -pe local 10
#$ -N ${sname}
#$ -hold_jid regMat-epimap-${HISTONE}-cut-${CUTOFF}-merge

echo "**** Job starts ****"
date

mkdir -p ${ROOTDIR}/${outdir}
mkdir -p ${ROOTDIR}/${outdir}/logs/

# Analyze results
module load R/3.3
Rscript ${ROOTDIR}/analyze_results_all.R -i "${HISTONE}" -t ${CUTOFF}

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
done
