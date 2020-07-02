#!/usr/bin/env Rscript

## load packages
library("GenomicFeatures")
library("AnnotationDbi")

## assign variables 
sourcegtf <- "/gpfs/data/kline-lab/ref/Gencode/gencode.v29.primary_assembly.annotation.gtf"
outfile <- "/gpfs/data/kline-lab/ref/Gencode/TxDb.gencode.v29.primary.sqlite"

# create TxDb object using GenomicFeatures package
TxDb.gencode.v29.primary <- makeTxDbFromGFF(sourcegtf, format="auto", dataSource="GENCODE", organism="Homo sapiens", taxonomyId=NA, circ_seqs=DEFAULT_CIRC_SEQS, chrominfo=NULL, miRBaseBuild=NA)

# save TxDb object for later use using AnnotationDbi package
saveDb(TxDb.gencode.v29.primary, outfile)

# exit when finished
quit()