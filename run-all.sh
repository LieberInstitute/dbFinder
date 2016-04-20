#!/bin/sh

## Usage
# sh run-all.sh shulha run8-v1.5.35 TRUE

# Define variables
EXPERIMENT=$1
PREFIX=$2
SKIP1=${3-"FALSE"}
SKIP6=${4-"FALSE"}
SKIP8=${5-"FALSE"}

mkdir -p ${EXPERIMENT}/CoverageInfo
mkdir -p ${EXPERIMENT}/derAnalysis
#mkdir -p ${EXPERIMENT}/regionMatrix
#mkdir -p ${EXPERIMENT}/regionMatrix-vs-DERs
#mkdir -p ${EXPERIMENT}/coverageToExon

if [[ $SKIP1 == "FALSE" ]]
then  
    sh step1-fullCoverage.sh ${EXPERIMENT}
fi
sh step2-makeModels.sh ${EXPERIMENT} ${PREFIX}
sh step3-analyzeChr.sh ${EXPERIMENT} ${PREFIX}
sh step4-mergeResults.sh ${EXPERIMENT} ${PREFIX}
sh step5-derfinderReport.sh ${EXPERIMENT} ${PREFIX}

#if [[ $SKIP6 == "FALSE" ]]
#then  
#    sh step6-regionMatrix.sh ${EXPERIMENT}
#fi
#sh step7-regMatVsDERs.sh ${EXPERIMENT} ${PREFIX}

#if [[ $SKIP8 == "FALSE" ]]
#then  
#    sh step8-coverageToExon.sh ${EXPERIMENT}
#fi

#sh step9-summaryInfo.sh ${EXPERIMENT} ${PREFIX}
