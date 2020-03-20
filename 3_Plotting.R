## P-value histogram

hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

## Nice single gene plot treatment vs control for gene with lowest p-value

plotCounts(dds, gene=which.min(res$padj), intgroup="treatment")

## Another single gene plot treatment vs control for gene with lowest p-value

d <- plotCounts(dds, gene=which.min(res$padj), intgroup="treatment",returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=treatment, y=count)) + geom_point(position=position_jitter(w=0.1,h=0)) + scale_y_log10(breaks=c(25,100,400))

## MA-plot: Points will be colored red if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.

plotMA(res05, ylim=c(-6,6))

## Plot with top DE gene circled:

plotMA(res, ylim = c(-5,5))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

##Identify:Use the function identify to interactively detect the row number of individual genes by clicking on the plot. One can then recover the gene identifiers by saving the resulting indices


idx <- identify(res05$baseMean, res05$log2FoldChange)
rownames(res05)[idx]

##Apeglm:Using apeglm adaptive t prior shrinkage estimator: only to show differentially expressed genes.

library(apeglm)
resApe <- lfcShrink(dds, coef=2, type="apeglm")
plotMA(resApe, ylim=c(-6,6))

##Ashr:Using ashr adaptive shrinkage estimator to fit a mixture of normal distributions to form the prior.

BiocManager::install("ashr")
#no

library("ashr")
resAsh <- lfcShrink(dds, coef=2, type="ashr")
plotMA(resAsh, ylim=c(-6,6))

##Comparison of shrinkage

resLFC <- lfcShrink(dds, coef="treatment_Treatment_vs_Control", type="apeglm")
resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")

##Three plots compare transformations

BiocManager::install("dplyr")
BiocManager::install("hexbin")
BiocManager::install("ggplot2")
#no

library("hexbin")
library("ggplot2")
library("dplyr")

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")  

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation) 

##Effects of transformations on the variance: Low count genes have a high standard deviation

ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

##Heatmap of the sample-to-sample distances:An overview over similarities and dissimilarities between samples

BiocManager::install("pheatmap")
library("pheatmap")
sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$treatment, vsd$samplename, sep="-")
colnames(sampleDistMatrix) <- paste(vsd$treatment, vsd$samplename, sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,col=colors)


