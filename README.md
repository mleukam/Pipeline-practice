# Pipeline Practice #

## Overview ##
This is a practice pipeline for sample WES data that was already analyzed by CRI. It is designed for two major tasks:
1. Preprocessing and alignment
2. Variant calling

The educational purposes of this pipeline are to:
1. Develop usable and portable code for preprocessing and variant calling
2. Practice using and troubleshooting the tools with real data
3. Gather information on the resources and time required for each analysis in Gardner HPC

The preprocessing pipeline follows very closely to GATK best practices, mostly using legacy tutorials from GATK 3.x and modifying where necessary for GATK 4.x

### Inputs ###
This pipeline assumes paired-end WES input in FASTQ format (file_1 and file_2). 

### Usage ###
Each step in the pipeline has a PBS script for batch submission to Gardner HPC at UChicago which points to a single BASH shell script. To use the scripts, ensure the scripts point to the correct directories and that all necessary scripts have had permissions changed to make executable:

`chmod +x <SCRIPTNAME>`

Go to the directory holding the executable PBS script and submit the job: 

`qsub <SCRIPTNAME.pbs>`

### Tool versions ###
* Java v1.8
* Picard tools version v2.8.1
* GATK v4.0.6.0
* BWA 0.7.17

## Reference sequences ##

### BWA ###
Reference alignment genome is version prepared by Heng Li (BWA author): GRCh38 with viral decoys and no alt contigs. Original filename: GRCh38_full_plus_hs38d1_analysis_set. See http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use for more info about the reference genome

### BaseRecalibrator ###
* dbSNP release 151: https://www.ncbi.nlm.nih.gov/projects/SNP/snp_summary.cgi
* Mills_and_1000G_gold_standard.indels.hg38.vcf from [GATK Resource Bundle](https://software.broadinstitute.org/gatk/download/bundle) (hg38)
* 1000G_phase1.snps.high_confidence.hg38.vcf from [GATK Resource Bundle](https://software.broadinstitute.org/gatk/download/bundle) (hg38)

### Reference indices ###
* GATK indices: 
	* GRCh38_full_plus_decoy.fa.fai (made with samtools faidx)
	* dbSNP, gold standard indels, and phase 1 VCF indices: `*.tbi` (made with IndexFeatureFile, see mlgatkindex.sh)
* Picard dictionary: GRCh38_full_plus_decoy.fa.dict (made with Picard tools)
* BWA indices (made with bwa index, see mlbwaindex.sh)

## Known issues ##
1. FASTQ to uBAM script (mlubam.sh) needs to have additional arguments for read group and platform information 
