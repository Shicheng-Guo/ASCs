for i in `ls *extendedFrags.fastq`
do
sed -n '1~4s/^@/>/p;2~4p'  $i > $i.fa
done
