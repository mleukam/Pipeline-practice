## Script to add missing information to a merged BAM after cleaning
## This should have been input when the uBAM was made, but wasn't 
## THis is required before base recalibration
## In the future this script should be unnecessary and is therefore a custom script
## Input is aligned BAM files after merger with uBAM (output of mlalign.sh)
## Output is marked BAM file with read information added
## Usage = ./mlfixbam.sh inside PBS script
## Built specifically for WES data that was transferred 8/2/18
## Designed for batch submission to Gardner HPC at UChicago
## First line = shebang to specify interpretor (bash)
## Before using, use chmod to make executable

#!/bin/bash

## Set script to fail if any command, variable, or output fails
set -euo pipefail

## Set IFS to split only on newline and tab
IFS=$'\n\t'

## Load compiler
module load java-jdk/1.8.0_92

## Load necessary modules
## NB: the path to picard.jar is automatically generated as an environment variable when picard tool module is loaded
## location of picard.jar = ${PICARD}
module load picard/2.8.1

## Navigate to directory with input aligned and merged BAM files
cd /scratch/mleukam/james_wes/

## Add missing information
java -Xmx8G -jar ${PICARD} AddOrReplaceReadGroups \
    INPUT=ABC-08-T_piped.bam \
    OUTPUT=ABC-08-T_addRG.bam \
    RGID=CCEMTANXX31 \
    RGLB=library1 \
    RGPL=illumina \
    RGPU=CCEMTANXX3ABC-08-T \
    RGSM=ABC-08-T \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true \
    TMP_DIR=/scratch/mleukam/temp

java -Xmx8G -jar ${PICARD} AddOrReplaceReadGroups \
    INPUT=ABC-09-T_piped.bam \
    OUTPUT=ABC-09-T_addRG.bam \
    RGID=CCEMTANXX31 \
    RGLB=library1 \
    RGPL=illumina \
    RGPU=CCEMTANXX3ABC-09-T \
    RGSM=ABC-09-T \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true \
    TMP_DIR=/scratch/mleukam/temp

java -Xmx8G -jar ${PICARD} AddOrReplaceReadGroups \
    INPUT=ABC-10-T_piped.bam \
    OUTPUT=ABC-10-T_addRG.bam \
    RGID=CCEMTANXX31 \
    RGLB=library1 \
    RGPL=illumina \
    RGPU=CCEMTANXX3ABC-10-T \
    RGSM=ABC-10-T \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true \
    TMP_DIR=/scratch/mleukam/temp

java -Xmx8G -jar ${PICARD} AddOrReplaceReadGroups \
    INPUT=ABC-11-T_piped.bam \
    OUTPUT=ABC-11-T_addRG.bam \
    RGID=CCEMTANXX41 \
    RGLB=library1 \
    RGPL=illumina \
    RGPU=CCEMTANXX4ABC-11-T \
    RGSM=ABC-11-T \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true \
    TMP_DIR=/scratch/mleukam/temp

java -Xmx8G -jar ${PICARD} AddOrReplaceReadGroups \
    INPUT=ABC-12-T_piped.bam \
    OUTPUT=ABC-12-T_addRG.bam \
    RGID=CCEMTANXX41 \
    RGLB=library1 \
    RGPL=illumina \
    RGPU=CCEMTANXX4ABC-12-T \
    RGSM=ABC-12-T \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true \
    TMP_DIR=/scratch/mleukam/temp

java -Xmx8G -jar ${PICARD} AddOrReplaceReadGroups \
    INPUT=AMP-06-T_piped.bam \
    OUTPUT=AMP-06-T_addRG.bam \
    RGID=CCEMTANXX41 \
    RGLB=library1 \
    RGPL=illumina \
    RGPU=CCEMTANXX4AMP-06-T \
    RGSM=AMP-06-T \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true \
    TMP_DIR=/scratch/mleukam/temp

java -Xmx8G -jar ${PICARD} AddOrReplaceReadGroups \
    INPUT=AMP-07-T_piped.bam \
    OUTPUT=AMP-07-T_addRG.bam \
    RGID=CCEMTANXX41 \
    RGLB=library1 \
    RGPL=illumina \
    RGPU=CCEMTANXX4AMP-07-T \
    RGSM=AMP-07-T \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true \
    TMP_DIR=/scratch/mleukam/temp

java -Xmx8G -jar ${PICARD} AddOrReplaceReadGroups \
    INPUT=AMP-08-T_piped.bam \
    OUTPUT=AMP-08-T_addRG.bam \
    RGID=CCEMTANXX51 \
    RGLB=library1 \
    RGPL=illumina \
    RGPU=CCEMTANXX5AMP-08-T \
    RGSM=AMP-08-T \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true \
    TMP_DIR=/scratch/mleukam/temp

java -Xmx8G -jar ${PICARD} AddOrReplaceReadGroups \
    INPUT=AMP-09-T_piped.bam \
    OUTPUT=AMP-09-T_addRG.bam \
    RGID=CCEMTANXX51 \
    RGLB=library1 \
    RGPL=illumina \
    RGPU=CCEMTANXX5AMP-09-T \
    RGSM=AMP-09-T \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true \
    TMP_DIR=/scratch/mleukam/temp

java -Xmx8G -jar ${PICARD} AddOrReplaceReadGroups \
    INPUT=AMP-10-T_piped.bam \
    OUTPUT=AMP-10-T_addRG.bam \
    RGID=CCEMTANXX51 \
    RGLB=library1 \
    RGPL=illumina \
    RGPU=CCEMTANXX5AMP-10-T \
    RGSM=AMP-10-T \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true \
    TMP_DIR=/scratch/mleukam/temp

java -Xmx8G -jar ${PICARD} AddOrReplaceReadGroups \
    INPUT=AMP-11-T_piped.bam \
    OUTPUT=AMP-11-T_addRG.bam \
    RGID=CCEMTANXX51 \
    RGLB=library1 \
    RGPL=illumina \
    RGPU=CCEMTANXX5AMP-11-T \
    RGSM=AMP-11-T \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true \
    TMP_DIR=/scratch/mleukam/temp

java -Xmx8G -jar ${PICARD} AddOrReplaceReadGroups \
    INPUT=AMP-12-T_piped.bam \
    OUTPUT=AMP-12-T_addRG.bam \
    RGID=CCEMTANXX61 \
    RGLB=library1 \
    RGPL=illumina \
    RGPU=CCEMTANXX6AMP-12-T \
    RGSM=AMP-12-T \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true \
    TMP_DIR=/scratch/mleukam/temp
