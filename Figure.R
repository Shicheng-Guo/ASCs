data<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/ASCs/master/extdata/vdjtools.basicstats.txt",head=T,sep="\t")
head(data)
idx1<-unlist(lapply(data[,1],function(x) substr(x,1,1)))
idx2<-unlist(lapply(colnames(input),function(x) unlist(strsplit(x,"[_]"))[2]))
idx<-paste(idx1,idx2,sep="_")
dim(data)
par(mfrow=c(3,4),mar=c(2,2,3,1))
for(i in 3:ncol(data)){
  boxplot(data[,i]~idx,col=2:5,main=colnames(data)[i])
}

data<-read.table("vdjtools.segments.wt.J.txt",head=T,sep="\t")
head(data)
idx1<-unlist(lapply(data[,1],function(x) substr(x,1,1)))
idx2<-unlist(lapply(colnames(input),function(x) unlist(strsplit(x,"[_]"))[2]))
idx<-paste(idx1,idx2,sep="_")
dim(data)
par(mfrow=c(2,4),mar=c(2,2,3,1))
for(i in 3:ncol(data)){
  boxplot(data[,i]~idx,col=2:5,main=colnames(data)[i])
}

data<-read.table("vdjtools.segments.wt.V.txt",head=T,sep="\t")
head(data)
dim(data)
idx1<-unlist(lapply(data[,1],function(x) substr(x,1,1)))
idx2<-unlist(lapply(colnames(input),function(x) unlist(strsplit(x,"[_]"))[2]))
idx<-paste(idx1,idx2,sep="_")
dim(data)
par(mfrow=c(4,4),mar=c(2,2,3,1))
for(i in 3:ncol(data)){
  if(sd(data[,i])>0.0055){
    boxplot(data[,i]~idx,col=2:5,main=colnames(data)[i],ylim=c(0,0.12))
  }
}


s1<-read.table("S2.txt",head=T,sep="\t")
s2<-read.table("vdjtools.m.txt",head=T,sep="\t")
s<-s2[match(s1[,1],gsub(".sed","",s2[,1])),]
write.table(s,file="vdjtools.m2.txt",col.names = T,row.names = F,sep="\t",quote=F)

vdjtools Convert -S mixcr -m vdjtools.m2.txt  metadata.txt
vdjtools CalcSegmentUsage -m metadata.txt vdjtools


data1<-read.table("vdjtools.spectratype.aa.wt.txt",head=T,sep="\t",check.names = F)
Len<-c()
for(i in 1:nrow(data1)){
  len<-sum(data1[i,7:22]*as.numeric(colnames(data1)[7:22]))
  Len<-c(Len,len)
}

tapply(Len,paste(data1$SampeType,data1$CellType,sep="_"),mean)
table(paste(data1$SampeType,data1$CellType,sep="_"))

pdf("Len.boxplot.pdf")
boxplot(Len~paste(data1$CellType,data1$SampeType,sep="_"),cex.axis=0.6,col=c(3,2,3,2,3,2,3,2))
dev.off()

x<-boxplot(Len~paste(data1$CellType,data1$SampeType,sep="_"),cex.axis=0.6,col=c(2,3,2,3,2,3,2,3))
colnames(x$stats)<-x$names
y<-x$stats
P<-c()
for(i in c(1,3,5,7)){
  p<-t.test(y[,i],y[,i+1])$p.value
  P<-c(P,p)
}
P
data2<-read.table("vdjtools.celltype.segments.wt.J.txt",head=T,sep="\t",check.names = F)
data2$IGHJ6
input<-data.frame(data1,Len,data2)

pdf("lenJH6.pdf")
plot(x=input$Len,y=input$IGHJ6,col=as.numeric(input$SampeType)+1,pch=as.numeric(input$CellType),cex=1.5)
dev.off()



data<-read.table("vdjtools.celltype.segments.wt.V.txt",head=T,sep="\t",check.names = F)
P<-c()
for(i in 6:ncol(data)){
  p<-t.test(data[,i]~data[,4])$p.value
  P<-c(P,p)
}
names(P)<-colnames(data)[6:ncol(data)]
qqplot(P)

data<-read.table("vdjtools.celltype.segments.wt.J.txt",head=T,sep="\t",check.names = F)
P<-c()
for(i in 6:ncol(data)){
  p<-t.test(data[,i]~data[,4])$p.value
  P<-c(P,p)
}
names(P)<-colnames(data)[6:ncol(data)]
qqplot(P)

meta<-read.table("metadata.txt",head=T,sep="\t")
for(i in unique(meta$SampleID)){
input<-subset(meta,SampleID==i)  
write.table(input,file=paste("metadata",i,"txt",sep="."),col.names =T,row.names = F,quote=F,sep="\t")
}

meta<-read.table("metadata.txt",head=T,sep="\t")
for(i in unique(meta$SampleType)){
  input<-subset(meta,SampleType==i)  
  write.table(input,file=paste("metadata",i,"txt",sep="."),col.names =T,row.names = F,quote=F,sep="\t")
}

library(reshape2)
library("ggplot2")
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/rheumatology/SLE/BCR/vdj")
Hs<-read.table("vdj.H.intersect.batch.aa.txt",head=T,sep="\t")
x<-subset(Hs,X1_SampleID==X2_SampleID &(X1_CellType=="BND" |X2_CellType=="BND" ))$R
Ss<-read.table("vdj.S.intersect.batch.aa.txt",head=T,sep="\t")
y<-subset(Ss,X1_SampleID==X2_SampleID &(X1_CellType=="BND" |X2_CellType=="BND" ))$R
input<-melt(data.frame(H=x,S=y))
input<-melt(data.frame(H=log(x,10),S=log(y,10)))
head(input)
wilcox.test(value~variable,input)
t.test(value~variable,input)
boxplot(value~variable,input)

pplot<-ggplot(input, aes(x=variable, y=value)) + 
  geom_boxplot(outlier.shape=NA,colour=c("green","red"),fill=c("green","red"))+
  geom_jitter(position=position_jitter(width=.1, height=0))+
  scale_y_continuous(name = "Log(D,10)")+
  scale_x_discrete(name = "") +
  theme_bw()
  add_pval(pplot, pairs = list(c(1, 2)), test='wilcox.test')
ggsave("BND.boxplot.ggplot.pdf")

