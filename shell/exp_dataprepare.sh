## load modules
module load gcc/6.2.0
module load R/3.4.1
module load samtools/1.6.0

## define variables
infile=$1
echo "input file is "${infile}

## navigate to directory containing perl scripts
cd /scratch/mleukam/EXCAVATOR2

## run excavator data prepare script
## formatted input list found in inputs folder in github repo - copied into shared lab space
perl EXCAVATORDataPrepare.pl ${infile} \
--processors 12 \
--target MyTarget_w50000 \
--assembly hg19
