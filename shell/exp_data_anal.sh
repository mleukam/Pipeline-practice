## Experimental Data Preparation for EXCAVATOR2 (version 2)

echo START `date`

## load modules

module load gcc/6.2.0 R/3.4.1 samtools/1.6.0 bedtools/2.26.0 perl/5.24.0

## define variables

project=dataAnalysis
threads=6
mode='pooling'
assembly=hg19
targetName=S04380110_V5_w50k

excavator2PATH=/gpfs/data/kline-lab/rbao_test/software/excavator2/1.1.2
projPath=/gpfs/data/kline-lab/rbao_test

targetBED=/gpfs/data/kline-lab/rbao_test/S04380110_V5_Covered.bed
expAnal=/gpfs/data/kline-lab/inputs/ExperimentalFileAnalysis.subset.w50K.txt

## must cd into excavator2PATH and run from within the installation dir...
cd $excavator2PATH

echo -e "project = $project\nassembly = $assembly\ntargetName = $targetName"

## loop through each sample per project
echo "[ " `date` " ] EXCAVATORDataAnalysis.pl"
echo -e "expAnal= $expAnal"

perl EXCAVATORDataAnalysis.pl $expAnal \
--processors ${threads} \
--target ${targetName} \
--assembly ${assembly} \
--output /scratch/mleukam/Results_Duke_w50K_2 \
--mode $mode


echo END `date`
