#!/usr/bin/env Rscript

# requires R version 3.5.0

# load packages
library("dplyr")
library("readr")
library("tximport")
library("GenomicFeatures")
library("AnnotationDbi")
library("edgeR")

# import input files
TxDb.gencode.v29.primary <- loadDb("/gpfs/data/kline-lab/ref/Gencode/TxDb.gencode.v29.primary.sqlite")
kallisto.files <- read_table("/gpfs/data/kline-lab/inputs/abundance_tsv.tsv", col_names=FALSE)
metadata <- read_csv("/gpfs/data/kline-lab/inputs/rna_meta_table.csv", col_names=FALSE)

# select columns to make tx2gene
k <- keys(TxDb.gencode.v29.primary, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(TxDb.gencode.v29.primary, k, "GENEID", "TXNAME")

# convert filenames and samples to vector
kallisto.files <- dplyr::pull(kallisto.files[,1])
sample.names <- dplyr::pull(metadata[,1])
names(kallisto.files) <- sample.names

# how to work around ignoreAfterBar
# https://www.biostars.org/p/322297/

# import reabdundance files for use with edger
txi.edger <- tximport(kallisto.files, 
	type = "kallisto", 
	tx2gene = tx2gene, 
	ignoreAfterBar = TRUE)

# import abundance files for use with limma
txi.limma <- tximport(kallisto.files,
	type = "kallisto",
	tx2gene = tx2gene,
	countsFromAbundance = "lengthScaledTPM",
	ignoreAfterBar = TRUE)

# save tximport objects for later use
saveRDS(txi.edger, "/gpfs/data/kline-lab/inputs/txi.edger.rds")
saveRDS(txi.limma, "/gpfs/data/kline-lab/inputs/txi.limma.rds")