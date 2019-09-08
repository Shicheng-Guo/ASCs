library("reshape2")
library("ggplot2")
library("ggpval")
############################################################################
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

############################################################################
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/rheumatology/SLE/BCR/vdj")
Hs<-read.table("vdj.H.intersect.batch.aa.txt",head=T,sep="\t")
x<-subset(Hs,X1_SampleID==X2_SampleID &(X1_CellType=="NAIVE" |X2_CellType=="NAIVE" ))$D
Ss<-read.table("vdj.S.intersect.batch.aa.txt",head=T,sep="\t")
y<-subset(Ss,X1_SampleID==X2_SampleID &(X1_CellType=="NAIVE" |X2_CellType=="NAIVE" ))$D
input<-melt(data.frame(H=x,S=y))
input
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
ggsave("Naive.logD.boxplot.ggplot.pdf")

Hs<-read.table("vdj.H.intersect.batch.aa.txt",head=T,sep="\t")
x<-subset(Hs,X1_SampleID==X2_SampleID &(X1_CellType=="NAIVE" |X2_CellType=="NAIVE" ))$R
Ss<-read.table("vdj.S.intersect.batch.aa.txt",head=T,sep="\t")
y<-subset(Ss,X1_SampleID==X2_SampleID &(X1_CellType=="NAIVE" |X2_CellType=="NAIVE" ))$R
input<-melt(data.frame(H=x,S=y))
head(input)
wilcox.test(value~variable,input)
t.test(value~variable,input)
boxplot(value~variable,input)
pplot<-ggplot(input, aes(x=variable, y=value)) + 
  geom_boxplot(outlier.shape=NA,colour=c("green","red"),fill=c("green","red"))+
  geom_jitter(position=position_jitter(width=.1, height=0))+
  scale_y_continuous(name = "Pearson correlation(R)")+
  scale_x_discrete(name = "") +
  theme_bw()
add_pval(pplot, pairs = list(c(1, 2)), test='wilcox.test')
ggsave("Naive.R.boxplot.ggplot.pdf")
