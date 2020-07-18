library(tximport)
library(tximportData)
library(DESeq2)
library(data.table)
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library(apeglm)
library("genefilter")

dir <- system.file("extdata", package = "tximportData")
list.files(dir)

samples <- read.table(file.choose())
samples$V1

# establish the folder for the abundance files and read them into a kallisto matrix
dir <- "/Users/mattsdwatson/mouse_de"
list.files(dir)
folders <- paste(samples$V1, "_kallisto", sep = "")
folders
files <- file.path(dir, folders, "abundance.tsv")
files
names(files) <- paste0(samples$V1)
kallisto_matrix <- tximport(files, type = "kallisto", ignoreAfterBar = TRUE, txOut = TRUE)
kallisto_matrix

sample_info <- data.frame(condition = c(rep("embryonic", 4), rep("immune", 2)))
row.names(sample_info) <- colnames(kallisto_matrix$counts)
sample_info

dds <- DESeqDataSetFromTximport(kallisto_matrix, sample_info, ~condition)
nrow(dds)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)

differential <- DESeq(dds)

summary(differential)

vsd <- vst(differential, blind = FALSE)
head(assay(vsd), 10)

# sample distance (euclidean) among samples
sampleDists <- dist(t(assay(vsd)))
sampleDists

sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# heatmap with the poisson distance among the samples
poisd <- PoissonDistance(t(counts(dds)))
pois_matrix <- as.matrix(poisd$dd)
rownames(pois_matrix) <- row.names(sample_info)
pheatmap(pois_matrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

# PCA for the different groups
plotPCA(vsd)

res <- results(differential)
res

res_shrink <- lfcShrink(differential, contrast=c("condition","embryonic","immune"), type = "apeglm")

ordered_p_value <- res[order(res$pvalue),]
summary(ordered_p_value)

# check for the number of adjusted p values
sum(res$padj < 0.05, na.rm=TRUE)

plotMA(res, ylim=c(min(res$log2FoldChange),max(res$log2FoldChange)), alpha = 0.05)
topgenes <- row.names(ordered_p_value[1:2,])
with(res[topgenes, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=1, lwd=2)
  text(baseMean, log2FoldChange, topgenes, pos=2, col="dodgerblue")
})

# view the top gene by significant differential expression

topGene <- rownames(res)[which.min(res$padj)]

#graph the count plot for the top 10 most differentially expressed genes
topgenes <- row.names(ordered_p_value[1:10,])
par(mfrow = c(3,4))
for (i in topgenes) {
  plotCounts(dds, gene = i, intgroup = kallisto_matrix$counts)
}

dev.off()
# histogram of p values
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "black", border = "white", main = "Histogram of p-values (DE)")

# gene clustering
# top 25 genes with the highest variance
# use the conditions as well as the sample names for comparison
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 25)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
pheatmap(mat, annotation_col = sample_info)

# volcano plot of log fold change and log p value

res$gene_significance <- ifelse(res$padj < 0.05, "significant", "not significant")
library(ggplot2)
data_frame <- as.data.frame(res)
qplot(log2FoldChange, -log(res$padj), data = data_frame, col = gene_significance)

# save the most variant genes to a csv for a Biomart search
genes_to_search <- c(row.names(mat))
write.csv(genes_to_search, "genes.csv")
