cd /gpfs/home/guosa/hpc/rheumatology/SLE/BCR/vdj
mkdir temp
for i in $(cat fq.txt)
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo mixcr align -f -s hsa -p rna-seq $i\_R1.fq $i\_R2.fq $i.vdjca -r $i  >>$i.job
echo mixcr assemble -f --write-alignments $i.vdjca $i.clna >>$i.job
echo mixcr assembleContigs -f $i.clna $i.clns >>$i.job
echo mixcr exportClones -f $i.clns $i.txt >>$i.job
qsub $i.job
done

cd /gpfs/home/guosa/hpc/rheumatology/SLE/BCR/vdj
mkdir temp
for i in $(cat fq.txt)
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo bbmerge.sh in1=$i\_R1.fq in2=$i\_R2.fq out=../imgt/$i.fq outu1=../imgt/$i.u1map outu2=../imgt/$i.u2map >>$i.job
qsub $i.job
done

for i in `ls *.fq`
do
sed -n '1~4s/^@/>/p;2~4p' $i > $i.fasta
echo $i
done

for i in `ls *.fasta`
do
seqtk sample -s seed=110 $i 150000 > $i.imgt
echo $i
done

for i in $(cat fq.txt)
do
sed -e 's/\s,\+/\t/g' -e 's/,\+\s/\t/g' $i.txt | sed -e 's/;,\+/;/g' > $i.sed
done

vdjtools Convert -S mixcr -m vdjtools.m.txt  metadata.txt
vdjtools CalcSegmentUsage -p -m metadata.txt vdjtools

vdjtools CalcSegmentUsage -f CellType -p -m metadata.txt vdjtools.celltype
vdjtools CalcSegmentUsage -f SampeName -p -m metadata.txt vdjtools.samplename
vdjtools CalcSegmentUsage -f SampeType -p -m metadata.txt vdjtools.sampletype

vdjtools CalcSpectratype -a -m metadata.txt vdjtools
vdjtools CalcPairwiseDistances -p -m metadata.txt vdjtools

for i in `ls metadata.*.txt`
do
vdjtools CalcPairwiseDistances -p -m $i vdj.$i
done
vdjtools CalcPairwiseDistances -p -m metadata.S.txt vdj.S
vdjtools ClusterSamples -p vdj.S  vdj.S
vdjtools CalcPairwiseDistances -p -m metadata.H.txt vdj.H
vdjtools ClusterSamples -p vdj.H  vdj.H

sed -e 's/\s,+/\t/g' -e 's/,+\s/\t/g' 16S1710LF01.txt | sed -e 's/;,+/;/g' >16S1710LF01.txt.sed.txt
sed -e 's/\s,\+/\t/g' -e 's/,\+\s/\t/g' 16S1710LF01.txt | sed -e 's/;,\+/;/g' > 16S1710LF01.txt.sed.txt

vdjtools PlotFancySpectratype [options] sample.txt output_prefix
vdjtools Convert -S mixcr -m vdjtools.txt  metadata.txt
qqplot(P)
vdjtools PlotFancyVJUsage vdjtools.segments.wt.V.txt 16S1710LF01
vdjtools CalcPairwiseDistances -p -m metadata.txt vdjtools

vdjtools Convert -S mixcr 16S1710LF01.txt.sed.txt  16S1710LF01.vdjtools


bbmerge.sh in1=<read1> in2=<read2> out=<merged reads> outu1=<unmerged1> outu2=<unmerged2>

for i in `ls *extendedFrags.fastq`
do
sed -n '1~4s/^@/>/p;2~4p'  $i > $i.fa
done


cp 1st/* vdj/
cp 2nd/* vdj/
gunzip *.gz 
http://ccb.jhu.edu/software/hisat2/dl/hisat2-2.1.0-source.zip
unzip hisat2-2.1.0-source.zip
cd /gpfs/home/guosa/hpc/tools/hisat2-2.1.0
make


setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/chol/16B1212A-2")
data= read_excel("methylation.xlsx",sheet = 2)
data= as.data.frame(data)
rowname<-apply(data.frame(data$Target,as.character(data$GenomePosition)),1,function(x) gsub(" ","",paste(x[1],x[2],sep="")))
data[1:12,1:12]
methdata<-data.matrix(data[,c(12:180)])
rownames(methdata)<-rowname
genesymbol= unlist(lapply(data$Target, function(x) strsplit(as.character(x),"_")[[1]][1]))
head(rownames(methdata))
head(colnames(methdata))

cat 180205LSJ32.fq.extendedFrags.fastq | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > $i.fa
