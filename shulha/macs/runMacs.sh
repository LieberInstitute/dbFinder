#!/bin/sh

## Usage
# sh runMacs.sh

# Directories
WDIR=/dcs01/ajaffe/Brain/derRuns/derChIP/shulha/macs

# Construct shell files
for sample in c1  c10  c11  c12  c13  c14  c15  c16  c17  c18  c19  c20  c21  c22  c23  c24  c25  c26  c27  c28  c29  c2p  c3  c30  c31  c32N  c33N  c34N  c35N  c3N  c4  c5  c6  c7  c8  c9
do
    sname="macs-${sample}"
    echo "Creating script ${sname}"
    cat > ${WDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=1G,h_vmem=3G,h_fsize=10G
#$ -N ${sname}

echo "**** Job starts ****"
date

## Create logs dir
mkdir -p ${WDIR}/logs

# merge results
cd ${WDIR}
module load macs/2.1.0
macs2 callpeak -t /dcs01/ajaffe/ChIPseq/Shulha2013/BED/${sample} -c /dcs01/ajaffe/ChIPseq/Shulha2013/BED/c28in --tsize 36 --bw 230 -n ${sample}_macs_out

# Move log files into the logs directory
mv ${WDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date
EOF
    call="qsub .${sname}.sh"
    echo $call
    $call
done
