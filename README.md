# Pipeline-practice

Practice pipeline for WES data from James' paper that was already analyzed by CRI. Should be easily modifiable for actual data analysis when primary data is ready.

Assumes paired-end WES input in FASTQ format (file_1 and file_2)
Each shell script has a partner PBS script to submit as a batch process to Gardner HPC at UChicago. To use the scripts, ensure the directories are correct that point to the shell scripts and that all necessary scripts have had permissions changed to make executable. Go to the directory holding the executable PBS script and submit the job using the "qsub" command.

Preprocessing pipeline follows very closely to GATK best practices, mostly using legacy tutorials from GATK 3.x and modifying where necessary for GATK 4.x

Picard tools version 2.8.1 is only version on Gardner HPC so the scripts are optimized for this version.

Most of the shell scripts use a design where all the file endings of interest in a particular directory are gathered into an array, and then the file endings are dropped to get individual sample names. The sample names are used to drive the "for" loop that iterates some task over the files of interest and to name the output. These will have to be carefully adjusted for future datasets to ensure that all filenames follow the same format.
