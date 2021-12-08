library(limma)
library(DESeq2)
library(sva)
library(dplyr)
library("RColorBrewer")

rm(list = ls())
setwd("D:/project/LDN-VC-endoderm/rna-seq/all_raw_count")
count<-read.table("RNA_oursampledata_count.txt",header = T,row.names = 1)
head(count)
count<-count[,-1:-4]
count<-count[,c(1,2,3,12,13,14,15,16,17,4,5,6,7,8,9,10,11,18,19)]
colnames(count)<-c("Day0-1","Day0-2","Day0-3","Day1-1","Day1-2","Day1-3",
                   "Vc-Day3-1","Vc-Day3-2","Vc-Day3-3","LDN-Day3-1","LDN-Day3-2",
                   "LDN-Day3-3","DE1-1","DE1-2","DE2-1","DE2-2","DE2-3","TE-1","TE-2")
data<-read.csv("GSE143371_COUNT_RNAseq_for_individual_mutants.csv")
data<-data %>% distinct(Gene_Name,.keep_all = T)
#data<- filter(data,!duplicated(data$Gene_Name))
rownames(data)<-data[,2]
data<-data[,c(3:5,9:10,17:18)]
head(data)
getwd()
head(count)
count<-count[,1:12]
count<-count[,c(1:6,10:12)]


str(count)
str(data)
gene<-rownames(count)
data<-data[gene,]

data2<-read.table("../Nissim Benvenisty/lastdata.txt",header = T)
head(data2)
data2<-data2%>% distinct(external_gene_name,.keep_all = T)
rownames(data2)<-data2$external_gene_name
head(data2)
str(data2)
data2<-data2[gene,]
head(rownames(data2))
head(rownames(count))
data<-data[,-1]
data2<-data2[,-1]
head(data)
last<-cbind(count,data)
last<-cbind(last,data2)
head(last)
str(last)
last$ME1<-as.integer(last$ME1)
last$ME2<-as.integer(last$ME2)
last$DE1<-as.integer(last$DE1)
last$DE2<-as.integer(last$DE2)
last$hESC1<-as.integer(last$hESC1)
last$hESC2<-as.integer(last$hESC2)
head(last)
#last<-apply(last,2,as.numeric)
library(DESeq2)
#treat <- factor(c(rep("VC",9),rep("LDN",3),rep("DE",5),rep("TE",2)))
#last<-last[,-1:-3]
head(last)
type<- factor(c(rep("Day0",3),rep("Day1",3),rep("DE",3),rep("Day0",2),rep("DE",2),rep("ME",2),
                rep("ME",2),rep("Day0",2),rep("DE",2)))
batch <- factor(c(rep("Current_study",9),rep("Yilmaz",6),rep("Haswell",6)))
colData <- data.frame(row.names=colnames(last),type=type,batch=batch)
colData

dim(colData)
str(last)

countData <- last[apply(last,1,sum) > 1 , ] 
head(countData)
which(is.na(countData))
countdata<-na.omit(countData)
dds<-DESeqDataSetFromMatrix(countdata,colData, formula(~batch+type)) 
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
# log转换后的结果
rld <- rlog(dds, blind=FALSE)
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]
#pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))
#hc <- hclust(t(rlogMat) )
?hclust,method="pearson"

library(ggplot2)
vsd <- vst(dds)
head(vsd)
pca_data <- plotPCA(vsd, intgroup=c("batch","type"), returnData=T, ntop=5000)
percentVar <- round(100 * attr(pca_data, "percentVar"))
ggplot(pca_data, aes(PC1, PC2, shape =batch,color=type)) +
  geom_point(size=3) +
  #xlim(-12, 12) +
  #ylim(-10, 10) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

####
library(limma)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), c(colData$batch))
pca <- plotPCA(vsd, intgroup=c("batch"), returnData=T, ntop=5000)
percentVar <- round(100 * attr(pca, "percentVar"))
head(pca)
p<-ggplot(pca, aes(PC1, PC2, shape =batch,color=type)) +
  geom_point(size=3) +
  #xlim(-12, 12) +
  #ylim(-10, 10) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

p+theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
head(assay(vsd))
pcacount<-assay(vsd)
head(pcacount)
write.table(pcacount,"RNAseq-PCAcount.txt")

#########get rm batch count matrix############
pcacount2<-pcacount
pcacount<-pcacount2
head(pcacount)
str(pcacount)
colnames(pcacount)<-c("Day0.1","Day0.2","Day0.3","Day1.1","Day1.2","Day1.3","LDN.Day3.1","LDN.Day3.2","LDN.Day3.3",
                      "Y_hESC2","Y_hESC3","Y_DE1","Y_DE2","Y_ME1","Y_ME2",
                      "H_ME1","H_ME2","H_hESC1","H_hESC2","H_DE1","H_DE2")

gene<-c("NANOG","FOXO1","POU5F1","PRDM14","SOX2","GSC","TBX6","MIXL1","TBXT","GATA4","TBX3","PITX1","HNF4A",
        "GATA3","SOX17","HNF1B","FOXA2","GATA6","FOXA1","HHEX","PIP4K2C","HIPK2","SPTBN1","DSC2","FNDC3B","KIT",
        "CADM1","CST1","CD177","CXCR4","AHNAK","RHOBTB3","FZD8","SHISA2","ZFP36L2","ERRFI1","SDC4","ATP6V1G1",
        "LIX1L","THBS1","CDH1","PPIC")
genetype<-c(rep("Pluripotency",5),rep("Mesendoderm",4),rep("Endoderm_TF",11),rep("Endoderm_pathway",22))
rowanno<-data.frame(type=factor(genetype))
rownames(rowanno)<-gene
colanno = data.frame(sampletype = factor(c(rep("hESC",7),rep("Mesendoderm",3),rep("Endoderm",7),rep("Mesoderm",4))))
rownames(colanno)<-colnames(pcacount)
data<-pcacount[gene,]
str(data)
data<-data[,c(1:3,10:11,18:19,4:9,12:13,20,21,14:17)]
rownames(colanno)<-colnames(data)

data=t(scale(t(data),scale = T,center = T))
library(pheatmap)
bk <- c(seq(-1,-0.1,by=0.01),seq(0,1,by=0.01))
#color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),
          #colorRampPalette(colors = c("white","red"))(length(bk)/2)),
#legend_breaks=seq(-8,8,2),
#breaks=bk
pheatmap(data,cluster_cols = F,cluster_rows =T ,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         cellwidth = 10, cellheight = 10,
         annotation_row = rowanno,
         annotation_col = colanno,
         show_rownames=T,show_colnames=T)

pheatmap(data,cluster_cols = F,cluster_rows = F,
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),
                   colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         #color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         breaks=bk,
         show_rownames=T,show_colnames=T)
####单独画####
Day1<-c("SMAD4","ZFP36L1","NODAL","BMPR1A","EOMES","LDB1","SSBP3","LHX1","APELA")
data<-pcacount2[Day1,]
str(data)
data<-data[,c(1:3,10:11,18:19,4:9,12:13,20,21,14:17)]
rownames(colanno)<-colnames(data)
data=t(scale(t(data),scale = T,center = T))
pheatmap(data,cluster_cols = F,cluster_rows =T ,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         cellwidth = 10, cellheight = 10,
         #annotation_row = rowanno,
         #annotation_col = colanno,
         show_rownames=T,show_colnames=T)
#######Sample DEG cluster tree#####
dds <- DESeq(dds)
H1vsDE<-results(dds,contrast =c("type","Day0","DE") )
H1vsME<-results(dds,contrast =c("type","Day0","ME") )
DEvsME<-results(dds,contrast =c("type","DE","ME") )
H1vsDay1<-results(dds,contrast =c("type","Day0","Day1") )
H1vsDay1_diff_gene <-subset(H1vsDay1, padj < 0.05 & abs(log2FoldChange)> 1) 
H1vsDE_diff_gene <-subset(H1vsDE, padj < 0.05 & abs(log2FoldChange)> 1) 
H1vsME_diff_gene <-subset(H1vsME, padj < 0.05 & abs(log2FoldChange)> 1) 
DEvsME_diff_gene <-subset(DEvsME, padj < 0.05 & abs(log2FoldChange)> 1) 
H1vsDE_diff_gene<-rownames(H1vsDE_diff_gene)
H1vsME_diff_gene<-rownames(H1vsME_diff_gene)
DEvsME_diff_gene<-rownames(DEvsME_diff_gene)
H1vsDay1_diff_gene<-rownames(H1vsDay1_diff_gene)
DEG<-c(H1vsDE_diff_gene,H1vsME_diff_gene,DEvsME_diff_gene,H1vsDay1_diff_gene)
DEG<-unique(DEG)

data<-pcacount[DEG,]
str(data)
data<-data[,c(1:3,10:11,18:19,4:9,12:13,20,21,14:17)]
data=t(scale(t(data),scale = T,center = T))
library(pheatmap)
bk <- c(seq(-1,-0.1,by=0.01),seq(0,1,by=0.01))
pheatmap(data,cluster_cols = T,cluster_rows =T ,
         clustering_method = "ward",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         show_rownames=F,show_colnames=T)

endoderm_develop<-c("GATA4","GATA4","KIF16B","KIF16B","KIF16B","SMAD2","SMAD4","GATA4","GATA4","NKX2-1","PAX9","BMP4","SOX17","ONECUT1","DKK1","MED12","NODAL","SMAD3","PELO","GDF3","NOTCH1","EXT1","ARC","BMPR1A","EOMES","HOXC11","WNT8A","APELA","APELA","APELA","TGFB1","HDAC1","EPB41L5","LAMC1","MIXL1","MIXL1","MIXL1","MIXL1","LHX1","HNF1B","NOG","BPTF","LHX1","SOX7","SMAD2","SOX7","DUSP4","DUSP4","SOX17","TBX20","DKK1","EOMES","EOMES","EOMES","EOMES","CTNNB1","DUSP5","DUSP1","DUSP2","MIXL1","LHX1","LHX1","NOG")
endoderm_develop<-unique(endoderm_develop)
data<-pcacount[which(rownames(pcacount)%in% endoderm_develop),]
str(data)
data<-data[,c(1:3,10:11,18:19,4:9,12:13,20,21,14:17)]
data=t(scale(t(data),scale = T,center = T))
library(pheatmap)
bk <- c(seq(-1,-0.1,by=0.01),seq(0,1,by=0.01))
pheatmap(data,cluster_cols = T,cluster_rows =T ,
         clustering_method = "ward",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         show_rownames=T,show_colnames=T)

#####merge rep#######
pcacount<-as.data.frame(pcacount)
pcacount$Day0    <-apply(pcacount[,1:3], 1, mean)
pcacount$Day1    <-apply(pcacount[,4:6], 1, mean)
pcacount$LDN_Day3<-apply(pcacount[,7:9], 1, mean)
pcacount$Y_hESC  <-apply(pcacount[,10:11], 1, mean)
pcacount$Y_DE    <-apply(pcacount[,12:13], 1, mean)
pcacount$Y_ME    <-apply(pcacount[,14:15], 1, mean)
pcacount$H_ME    <-apply(pcacount[,16:17], 1, mean)
pcacount$H_hESC  <-apply(pcacount[,18:19], 1, mean)
pcacount$H_DE    <-apply(pcacount[,20:21], 1, mean)
str(pcacount)
count<-pcacount[,22:30]
str(count)
count<-count[,c(1,3,8,2,5,4,9,6,7)]
DEGcount<-count[DEG,]
library(psych)
pdf("DEGcount-cor.pdf")
pairs.panels(DEGcount,cex.cor=1.5)
dev.off()


library(pheatmap)
pcacount<-read.table("D:/project/LDN-VC-endoderm/figure and MS/RNAseq-PCAcount.txt")
colnames(pcacount)<-c("Day0.1","Day0.2","Day0.3","Day1.1","Day1.2","Day1.3","LDN.Day3.1","LDN.Day3.2","LDN.Day3.3",
                      "Y_hESC2","Y_hESC3","Y_DE1","Y_DE2","Y_ME1","Y_ME2",
                      "H_ME1","H_ME2","H_hESC1","H_hESC2","H_DE1","H_DE2")

liver<-c("HNF4A","HHEX","ALB","HNF6","HNF1A","UGT1A1","ASGR","AAT","CYP3A7","CYP3A4","CYP1A1")
lung<-c("P63","CK5","NGFR","ACTTUB","SFTPC","SFTPB","ABCA3","LAMP3","LPCAT1","SCGB1A1","PDPN","AGER","CAV1")
pancreas<-c("INSULIN","NKX2-2","NKX6.1","NGN3","PDX1","MAFA","NEUROD1","CHGA","SST","ABCC8","G6PC2","PCSK1","PCSK2")
data<-pcacount[which(rownames(pcacount)%in% liver),]
str(data)

data<-data[,c(1:3,10:11,18:19,4:9,12:13,20,21,14:17)]
data=t(scale(t(data),scale = T,center = T))
library(pheatmap)
bk <- c(seq(-1,-0.1,by=0.01),seq(0,1,by=0.01))
pheatmap(data,cluster_cols = T,cluster_rows =T ,
         clustering_method = "ward.D2",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         show_rownames=T,show_colnames=T)

