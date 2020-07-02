## Experimental Data Preparation for EXCAVATOR2 (version 2)

echo START `date`

## load modules

module load gcc/6.2.0 R/3.4.1 samtools/1.6.0 bedtools/2.26.0 perl/5.24.0

## define variables

project=redo_list_1
threads=12
window=50000
assembly=hg19
targetName=S04380110_V5_w50k

excavator2PATH=/gpfs/data/kline-lab/rbao_test/software/excavator2/1.1.2
projPath=/gpfs/data/kline-lab/rbao_test

targetBED=/gpfs/data/kline-lab/rbao_test/S04380110_V5_Covered.bed
source=$excavator2PATH/data/targets/${assembly}/excavator2-${assembly}-SourceTarget.txt
expPrepapre=/gpfs/data/kline-lab/inputs/ExperimentalFilePrepare.w50000.missing.txt

## must cd into excavator2PATH and run from within the installation dir...
cd $excavator2PATH

echo -e "project = $project\nassembly = $assembly\nagilent = $agilent\ntargetName = $targetName"

## loop through each sample per project
echo "[ " `date` " ] Parsing BAM and normalization: EXCAVATORDataPrepare.pl"
echo -e "expPrepapre = $expPrepapre"
perl EXCAVATORDataPrepare.pl $expPrepapre --processors $threads --target $targetName --assembly $assembly

echo END `date`
