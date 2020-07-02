## load modules
module load gcc/6.2.0
module load R/3.4.1
module load samtools/1.6.0

## navigate to directory containing perl scripts
cd /scratch/mleukam/EXCAVATOR2

## run excavator data prepare script
## formatted input list found in inputs folder in github repo - copied into shared lab space
perl EXCAVATORDataPrepare.pl /gpfs/data/kline-lab/inputs/ControlFilePrepare.w50000.txt \
--processors 12 \
--target MyTarget_w50000 \
--assembly hg19
