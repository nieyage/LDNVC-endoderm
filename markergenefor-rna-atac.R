getwd()
library(pheatmap)
library("RColorBrewer")
library(DESeq2)
library(gplots)
setwd("D:/atac/all_raw_count")
options(stringsAsFactors = FALSE)
rm(list = ls())
count=read.table("markeratac.txt",header = T)
head(count)
count<-count[,c(1,2,3,12,13,14,15,4,5,16,17,6,7,8,9,10,11)]
rownames(count)<-count[,1]
count<-count[,-1]
#batch<-factor(c(rep("a",2),rep("b",10),rep("c",4)))
count<-count[,-1:-8]
bat<-factor(c(rep("a",4),rep("b",4)))
condition<-factor(c(rep("Vc",2),rep("LDN",2),rep("DE",2),rep("FOXA2KO",2)))
colData <- data.frame(row.names=colnames(count),bat,condition)
colData
design<-~bat+condition
dds<-DESeqDataSetFromMatrix(count,colData,~bat) 
dds <- estimateSizeFactors(dds)
dds<- estimateDispersions(dds)
dds<- DESeq(dds)
norm_count<-counts(dds,normalized=TRUE)
######sample map####
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

#############
####对counts log 转换
rld <-rlogTransformation(dds)
#######方差稳定转换
vsd <- varianceStabilizingTransformation(dds)
rlogMat <- assay(rld)
vstMat <- assay(vsd)
###
normalizedCounts1 <- t( t(counts(dds)) / sizeFactors(dds) )
head(normalizedCounts1)
normalizedCounts2 <- counts(dds, normalized=T)
######相当于TPM
head(normalizedCounts2)
exprMatrix_rpm=as.data.frame(normalizedCounts2) 
head(exprMatrix_rpm)
exprMatrix_rpm=t(scale(t(exprMatrix_rpm),scale = T))
####limma remove batch effect####
library(ggplot2)
library("limma")
assay(rld) <- limma::removeBatchEffect(assay(rld), vsd$batch)
plotPCA(vsd, "batch")
pcadata<-plotPCA(rld, intgroup=c("batch"),returnData=TRUE)
ggplot(pcadata, aes(PC1, PC2 ,color=batch)) +geom_point(size=3) +
  coord_fixed()
ex_b_limma <- removeBatchEffect(exprMatrix_rpm,
                                batch = batch)
dim(ex_b_limma) 
head(ex_b_limma)

#######before rm #####
pcaData <- plotPCA(vsd, intgroup=c("bat"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=bat)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text(aes(label=name),vjust=2)
#######after rm #####
assay(vsd) <-limma::removeBatchEffect(assay(vsd), vsd$bat)
pcaData <- plotPCA(vsd, intgroup=c("bat"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=bat)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text(aes(label=name),vjust=2)

######data for heatmap#####
pheatmap(ex_b_limma, cluster_rows=TRUE, show_rownames=T,
         color = colorRampPalette(colors = c("blue","white","red"))(100),
         cluster_cols=T,cutree_rows =1,main = "ATAC marker gene Chromatin accessibility",
         show_colnames = T,legend_labels = "lll")
#############RNAseq##################

rm(list = ls())
setwd("D:/atac/rna-seq/all_raw_count")
count=read.table("RNA_DEdata_count.txt",header = T)
head(count)
rownames(count)<-count[,1]
count<-count[,-1:-2]
count<-count[,c(1,2,3,9,10,11,4,5,6,7,8)]
gene=read.table("markergene.txt",header = T)
gene<-gene$gene
count<-count[which(gene%in%rownames(count)),]
######DESEQ #########
batch<-factor(c(rep("a",6),rep("b",5)))
colData <- data.frame(row.names=colnames(count),batch)
colData
dds<-DESeqDataSetFromMatrix(count,colData,~batch) 
dds <- estimateSizeFactors(dds)
dds<- estimateDispersions(dds)
dds<- DESeq(dds)
norm_count<-counts(dds,normalized=TRUE)
rld <-rlogTransformation(dds)
#######方差稳定转换
vsd <- varianceStabilizingTransformation(dds)
rlogMat <- assay(rld)
vstMat <- assay(vsd)
####limma remove batch effect####
library(ggplot2)
library("limma")
exprMatrix_rpm=as.data.frame(norm_count) 
head(exprMatrix_rpm)
exprMatrix_rpm=t(scale(t(exprMatrix_rpm),scale = T))
ex_b_limma <- removeBatchEffect(exprMatrix_rpm,
                                batch = batch)
dim(ex_b_limma) 
head(ex_b_limma)

#######before rm #####
pcaData <- plotPCA(vsd, intgroup=c("batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=batch)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text(aes(label=name),vjust=2)
#######after rm #####
assay(rld) <- limma::removeBatchEffect(assay(rld), rld$batch)

assay(vsd) <-limma::removeBatchEffect(assay(vsd), vsd$bat)
pcaData <- plotPCA(rld, intgroup=c("batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=batch)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text(aes(label=name),vjust=2)

######data for heatmap#####
pheatmap(ex_b_limma, cluster_rows=TRUE, show_rownames=T,
         color = colorRampPalette(colors = c("blue","white","red"))(100),
         cluster_cols=T,cutree_rows =1,main = "RNA marker gene expression",
         show_colnames = T,legend_labels = "lll")
exprMatrix_rpm
pheatmap(exprMatrix_rpm, cluster_rows=TRUE, show_rownames=T,
         color = colorRampPalette(colors = c("blue","white","red"))(100),
         cluster_cols=T,cutree_rows =1,main = "RNA marker gene expression",
         show_colnames = T,legend_labels = "lll")
