## A20 Mutation Calls Version 2.0

##---------------------------------------

## README ##

## This is not meant to be a production pipeline script
## This was run in an interactive terminal on a local machine

# align short (50bp) *single end* whole genome reads from mouse tumor
# in addition to aligning, we will
# -- mark illumina adapters
# -- filter unmapped reads and multimapped reads
# -- index results of alignment and sort by coordinate
# -- merge results with unaligned bam to restore read group info and correct hard clipping
# -- remove duplicates after alignment

## pre-processing best practices
## https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery

## reference genome for this project
## use the same reference genome as the mouse genome project: ftp://ftp-mouse.sanger.ac.uk/ref/
## GRCm38_68

## installation of docker following instructions here:
## https://gatk.broadinstitute.org/hc/en-us/articles/360035889991
## installed docker from here:
## https://hub.docker.com/editions/community/docker-ce-desktop-mac/

## bioinfo conda environment from: https://www.biostarhandbook.com/install-software.html

###############
## DATA PREP ##
###############


## test installation
docker --version
## Docker version 19.03.8, build afacb8b

## downloaded gatk docker version from command line
docker pull broadinstitute/gatk:4.1.3.0

## test start up gatk container
docker run -it broadinstitute/gatk:4.1.3.0
gatk --help
exit

## unzip reference file
cd ~/projects/mouse
gunzip GRCm38_68.fa.gz

## before launching docker, index reference genome with bwa
## -a bwtsw specifies that we want to use the indexing algorithm that is capable of handling the whole human genome.
## load bioconda environment with bwa installed
## stay in working folder ~/projects/mouse
conda activate bioinfo

## index reference file for alignment
bwa index -a bwtsw GRCm38_68.fa

# return to base environment
conda activate

## mount a GATK docker instance with external folder mounted to see the data
docker run -v ~/projects/mouse:/gatk/mouse_data -it broadinstitute/gatk:4.1.3.0

## use Picard tools to create sequence dictionary
gatk --java-options "-Xmx16G" CreateSequenceDictionary \
-R mouse_data/GRCm38_68.fa 

## return to base environment
exit

## enter bioinfo environment
conda activate bioinfo

## check for fastqc installation
fastqc --help

## prepare sequence files and get quality report
## stay in working folder ~/projects/mouse
## unzip file and rename
gunzip -c JK-ST-DNA-A20_S1_L006_R1_001.fastq.gz > a20.fq

## make fastq folder and run fastqc
mkdir fastqc
fastqc -o fastqc a20.fq

# return to base environment
conda activate

## Convert fastq to uBAM
## Inputs: fq file from arguments
## Output: unaligned BAM
## NB: include necessary read information - this will need to be added by hand!
docker run -v ~/projects/mouse:/gatk/mouse_data -it broadinstitute/gatk:4.1.3.0

gatk --java-options "-Xmx16G" FastqToSam \
--FASTQ mouse_data/a20.fq \
--OUTPUT mouse_data/a20_unaligned.bam \
-SM A20 \
-LB coreWGS \
-PL illumina \
-PU K00242 \
-RG HW5FWBBXX.6 \
-DT 2018-09-30T00:00:00-0400 \
-CN UCHICAGOCORE

## Mark Illumina adapters
## Input: unaligned BAM
## Output: marked Illumina adapters
gatk --java-options "-Xmx16G" MarkIlluminaAdapters \
-I mouse_data/a20_unaligned.bam \
-O mouse_data/a20_markilluminaadapters.bam \
-M mouse_data/a20_markilluminaadapters_metrics.txt

## Save uBAM for later merger

#################################
## ALIGN EXPERIMENTAL SEQUENCE ##
#################################

## Initial alignment, merge with uBAM
## Input: marked Illumina adapters BAM
## Output: fastq
gatk --java-options "-Xmx16G" SamToFastq \
-I mouse_data/a20_markilluminaadapters.bam \
-FASTQ mouse_data/a20_reconstituted.fq \
-CLIPPING_ATTRIBUTE XT \
-CLIPPING_ACTION 2 \
-INTERLEAVE false \
-NON_PF true


# FileTruncatedException: Premature end of file: /gatk/mouse_data/a20_markilluminaadapters.bam

## get out of docker environment
exit

## load bioinfo conda environment
conda activate bioinfo

## align
## Input: fastq
## Output: aligned bam
bwa mem -t 4 -M -v 3 \
~/projects/mouse/GRCm38_68.fa \
~/projects/mouse/a20_reconstituted.fq >\
~/projects/mouse/a20_aligned.bam

## return to base environment
conda activate

## activate docker environment
docker run -v ~/projects/mouse:/gatk/mouse_data -it broadinstitute/gatk:4.1.3.0

## merge aligned bam with unaligned bam
## Input: aligned and unaligned bam
## Output: aligned bam with marked adapters, quality, and metadata preserved
gatk --java-options "-Xmx16G" MergeBamAlignment \
--ALIGNED_BAM mouse_data/a20_aligned.bam \
--UNMAPPED_BAM mouse_data/a20_unaligned.bam \
--OUTPUT a20_merged.bam \
-R mouse_data/GRCm38_68.fa \
--CREATE_INDEX true \
--ADD_MATE_CIGAR true \
--CLIP_ADAPTERS true \
--CLIP_OVERLAPPING_READS true \
--INCLUDE_SECONDARY_ALIGNMENTS true \
--MAX_INSERTIONS_OR_DELETIONS -1 \
--PRIMARY_ALIGNMENT_STRATEGY MostDistant \
--ATTRIBUTES_TO_RETAIN XS

## mark duplicates
gatk --java-options "-Xmx16G" MarkDuplicates \
--INPUT mouse_data/a20_merged.bam \
-OUTPUT mouse_data/a20_markduplicates.bam \
--METRICS_FILE mouse_data/a20_markduplicates_metrics.txt \
--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
--CREATE_INDEX true

## leave docker environment
exit

## activate base environement and navigate to mouse folder
conda activate 
cd ~/projects/mouse/

## Clean up: delete unused intermediate files
rm a20_unaligned.bam
rm a20_markilluminaadapters.bam
rm a20_reconstituted.bam
rm a20_merged.bam

## known site polymorphisms from balb/c mouse genome project
gunzip -c BALB_cJ.mgp.v5.snps.dbSNP142.vcf.gz > knownsites.snps.balbc.vcf
gunzip -c BALB_cJ.mgp.v5.indels.dbSNP142.normed.vcf.gz > knownsites.indels.balbc.vcf

# Clean up: delete intermediate files
rm BALB_cJ.mgp.v5.snps.dbSNP142.vcf.gz.tbi
rm BALB_cJ.mgp.v5.snps.dbSNP142.vcf.gz
rm BALB_cJ.mgp.v5.indels.dbSNP142.normed.vcf.gz.tbi
rm BALB_cJ.mgp.v5.indels.dbSNP142.normed.vcf.gz

## return to docker environment
docker run -v ~/projects/mouse:/gatk/mouse_data -it broadinstitute/gatk:4.1.3.0

## Base quality score re-alignment
## Input: marked duplicates BAM
## Output: BQSR table and recalibrated BAM (bqsr.bam)
gatk --java-options "-Xmx16G" BaseRecalibrator \
-R mouse_data/GRCm38_68.fa \
-I a20_markduplicates.bam \
--known-sites knownsites.snps.balbc.vcf \
--known-sites knownsites.indels.balbc.vcf \
-O a20_bqsr.table

gatk --java-options "-Xmx16G" ApplyBQSR \
-R mouse_data/GRCm38_68.fa \
-I a20_markduplicates.bam \
--bqsr-recal-file a20_bqsr.table \
-O a20_bqsr.bam

# leave docker environment
exit

#############################
## ALIGN GERMLINE SEQUENCE ##
#############################

## move downloaded BAM to working folder
## Balbc/J sequence downloaded from sanger mouse genome project:
## tp://ftp-mouse.sanger.ac.uk/current_bams

## return to docker environment
docker run -v ~/projects/mouse:/gatk/mouse_data -it broadinstitute/gatk:4.1.3.0

# convert bam to unaligned bam file
# include necessary read information - this will need to be added by hand!
gatk --java-options "-Xmx16G" RevertSam \
-I BALB_cJ.bam \
-O BALB_cJ_unaligned.bam \
--SANITIZE true \
--MAX_DISCARD_FRACTION 0.010 \
--ATTRIBUTE_TO_CLEAR XT \
--ATTRIBUTE_TO_CLEAR XN \
--ATTRIBUTE_TO_CLEAR AS \
--ATTRIBUTE_TO_CLEAR OC \
--ATTRIBUTE_TO_CLEAR OP \
--SORT_ORDER queryname \
--RESTORE_ORIGINAL_QUALITIES true \
--REMOVE_DUPLICATE_INFORMATION true \
--REMOVE_ALIGNMENT_INFORMATION true





#######################
## IDENTIFY VARIANTS ##
#######################



## generate VCF and bamout with mutect2
## can use four threads
java -Xmx16G -jar ${GATK} Mutect2 \
-R ${ref} \
-I ${tumor} \
-I ${normal} \
--tumor-sample ${sample} \
--normal-sample ${norm_name} \
--germline-resource ${knownsites2} \
--output ${sample}.vcf.gz \
-bamout tumor_normal_${sample}.bam \
--TMP_DIR ${tmpdir}
