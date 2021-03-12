#####QC Table#####
library("ggplot2")
setwd("D:/atac/ATACresults")
options(stringsAsFactors = FALSE)
rm(list = ls())
qc=read.csv("atacseq-data-information.csv",header = T)
head(qc)
qc
qc<-qc[-13,-17:-19]
qc<-qc[c(1,2,)]
rownames(qc)<-qc[,1]
names(qc)[1]<-"sample"
sample<-qc[,1]
sample<-sample[-13]
sample
####total reads#####
total<-qc[,2]
total
total<-data.frame(x=sample,y=total)
#total<-as.numeric(total)
head(total)
#class(total)
sample<-factor(sample,levels = c("D0-1", "D0-2" , "Vc-D1-1","Vc-D1-2", "LDN-D2-1", "LDN-D2-2" ,"Vc-D2-1", 
                      "Vc-D2-2","LDN-D3-1" ,"LDN-D3-2", "Vc-D3-1" , "Vc-D3-2"),ordered = T)

ggplot(data=total, mapping=aes(x=sample,y=total$y))+
  geom_bar(mapping = NULL,stat= 'identity', width=0.7,fill = 'steelblue')+
  labs(x="",y="total reads number",title="total reads")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_abline(intercept = 0, slope = 0,size=1,colour='black')
  )
#####mapped reads####
map<-qc[,16]
map=c("3202799","3338003","5273144","5267182","5377880","5362157","5428627","5554481","5937167","5955261","5505334","5386569")
map<-data.frame(x=sample,y=map)
ggplot(data=map, mapping=aes(x=sample,y=map$y))+
  geom_bar(mapping = NULL,stat= 'identity', width=0.7,fill = 'steelblue')+
  labs(x="",y="reads number",title="Number of Reads in blacklist regions")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_abline(intercept = 0, slope = 0,size=1,colour='black')
)

######RNA seq QC table######

setwd("D:/atac/rna-seq")
options(stringsAsFactors = FALSE)
rm(list = ls())
qc=read.csv("endoderm RNA-seq information.csv",header = T);
head(qc)
rownames(qc)<-qc[,1]
names(qc)[1]<-"sample"
sample<-qc[,1]
sample
sample<-factor(sample,levels = c("D0-1","D0-2","D0-3","Vc-D1-1",  "Vc-D1-2" , "Vc-D1-3" ,
                                 "LDN-D2-1" ,"LDN-D2-2", "LDN-D2-3" ,
                                 "Vc-D2-1" , "Vc-D2-2" , "Vc-D2-3",
                                 "LDN-D3-1", "LDN-D3-2" ,"LDN-D3-3" ,
                                 "Vc-D3-1" , "Vc-D3-2", "Vc-D3-3" ),ordered = T)

num<-qc[,2]
num<-data.frame(x=sample,y=qc$GC)
num
ggplot(data=num, mapping=aes(x=sample,y=num$y))+
  geom_bar(mapping = NULL,stat= 'identity', width=0.7,fill = 'steelblue')+
  labs(x="",y="GC rate",title="GC content")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_abline(intercept = 0, slope = 0,size=1,colour='black')
)
