## Install required packages:
BiocManager::install("DESeq2")
BiocManager::install("tximport")
BiocManager::install("tximportData")
BiocManager::install("AnnotationDbi")
BiocManager::install("AnnotationHub")
BiocManager::install("GenomicFeatures")
BiocManager::install("GenomicAlignments")
BiocManager::install("GenomicRanges")
BiocManager::install("GenomeInfoDb")


BiocManager::install("digest")
BiocManager::install("readr")
BiocManager::install("apeglm")
BiocManager::install("ashr")
BiocManager::install("truncnorm")
BiocManager::install("foreign")
BiocManager::install("lattice")
BiocManager::install("MASS")
BiocManager::install("Matrix")
BiocManager::install("mgcv")
BiocManager::install("survival")
BiocManager::install("hexbin")
BiocManager::install("dplyr")
BiocManager::install("ggplot2")
BiocManager::install("pheatmap")

##Load required libraries:
  
library("DESeq2")
library("tximport")
library("tximportData")
library("AnnotationDbi")
library("AnnotationHub")
library("GenomicFeatures")
library("GenomicAlignments")
library("GenomicRanges")
library("GenomeInfoDb")
library("AnnotationHub")
library("digest")
library("readr")

##Reading data

dir <- "/Users/<your_username>/desktop/RNA-seq_Differential_Expression/Demo_Data" 
samples <- read.table(file.path(dir, "Bias_samples.txt"), header = TRUE)
files <- file.path(dir, samples$directory, "quant.sf")
names(files) <- samples$samplename

##Read in a table that links transcripts to genes for this dataset and match gene Id with transcriptome

TxDb <- makeTxDbFromGFF(file = "/Users/<your_username>/desktop/RNA-seq_Differential_Expression/Demo_Data/gencode.vM16.annotation.gtf.gz")
k <- keys(TxDb, keytype = "TXNAME")
tx2gene <- select(TxDb, k, "GENEID", "TXNAME")
head(tx2gene)

##Import the necessary quantification data for DESeq2 using the tximport function

txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)

##construct a DESeqDataSet from the txi object and sample information in samples

library("DESeq2")
ddsTxi <- DESeqDataSetFromTximport(txi.salmon,
                                   colData = samples,
                                   design = ~ treatment)

##Sets "Control" as reference level (H0)

ddsTxi$treatment <- relevel(ddsTxi$treatment , ref = "Control")


## Pre-filtering options: pre-filtering to keep only rows that have at least 10 reads total
keep <- rowSums(counts(ddsTxi)) >= 10
dds <- ddsTxi[keep,]


##. DESeq function
dds <- DESeq(dds)
res <- results(dds)
res

## Sort by p-value, smallest to largest

resOrdered <- res[order(res$pvalue),]

## Count How many adjusted p-values were less than 0.05


sum(res$padj < 0.05, na.rm=TRUE)

## Run FDR at alpha=0.05. By default the argument alpha is set to 0.1. If the adjusted p value cutoff will be a value other than 0.1, alpha should be set to that value(e.g. 0.05)

res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

## Simple summary as to how many genes are significantly DE after FDR:
  
res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)

## If we want to raise the log2 fold change threshold, so that we test for genes that show more substantial changes due to treatment, we simply supply a value on the log2 scale. For example, by specifying lfcThreshold = 1, we test for genes that show significant effects of treatment on gene counts more than doubling or less than halving.

resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)






