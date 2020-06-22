# List of useful commands for shell/gardner

# Login
ssh mleukam@gardner.cri.uchicago.edu

# obtain a unique interactive node
qsub -I

# logout of unique interactive node
logout

# show available modules
module avail

# load module (example is bowtie)
# if no version is specified, loads default (marked with D on listing)
module load bowtie

# unzip .gz files and move to new directory (use when you are in the directory of the source file)
# remember to change target directory to intended file name
gunzip -c ABC-09-T_1.fq.gz > /scratch/mleukam/james_wes/ABC-09-T_1.fq

# run fastqc and output to a specific directory
# fastqc [-o (specify output)] [output directory] [input fastq file]
# note: fastqc won't create the directory so that has to be done beforehand
# creates zip file and html file in the target directory
fastqc -o fastqc/ABC-06-T_2_fastqc ABC-06-T_2.fq

# Good DNA alignment tools:
# bwa or bowtie

# Good converstion tool of BAM --> fastq
# SAMtools

## GARDNER TIPS:

# For batch jobs 
## = comments
#PBS is directive to scheduler

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#PBS -N <job_name>

## define the shell

#PBS -S /bin/bash

## request resources, single node, single processer, RAM
#PBS -l walltime=8:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=4gb

## specify the standard output file (rather than printing to screen)
#PBS -o /home/mleukam/output.log

## specify the standard error (rather than printing to screen)
#PBS -e /home/mleukam/error.log

## load compiler
module load gcc/6.2.0

## load programs
module load R/3.4.1

## executable section



## when you submit a job, it will give you a job number. Can check the status of the job with qstat (queue status)

qstat
watch qstat

## to show how much time is remaining
showq -u $USER
