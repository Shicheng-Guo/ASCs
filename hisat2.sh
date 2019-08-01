cp 1st/* vdj/
cp 2nd/* vdj/
gunzip *.gz 
http://ccb.jhu.edu/software/hisat2/dl/hisat2-2.1.0-source.zip
unzip hisat2-2.1.0-source.zip
cd /gpfs/home/guosa/hpc/tools/hisat2-2.1.0
make

mkdir indexes
cd indexes
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38_tran.tar.gz

