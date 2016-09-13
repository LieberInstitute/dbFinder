#!/bin/bash
#$ -cwd
#$ -m e
#$ -N makeBai-USC
#$ -pe local 10

## Create logs dir
mkdir -p logs

echo "**** Job starts ****"
date

module load samtools/1.1

for i in /dcl01/lieber/ajaffe/psychENCODE_Data/USC_U01MH103346/BAM/*.bam
    do echo $i
    ## Figure out number of mapped reads
    samtools view -c -F 4 $i
done

## Create .bam.bai files
find /dcl01/lieber/ajaffe/psychENCODE_Data/USC_U01MH103346/BAM -name '*.bam' | parallel --jobs 10 samtools index

echo "**** Job ends ****"
date

## Move log files
mv makeBai-USC.* logs/
