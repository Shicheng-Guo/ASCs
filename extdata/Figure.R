
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
