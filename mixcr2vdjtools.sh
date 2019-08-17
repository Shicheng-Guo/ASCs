for i in $(cat fq.txt)
do
sed -e 's/\s,\+/\t/g' -e 's/,\+\s/\t/g' $i.txt | sed -e 's/;,\+/;/g' > $i.sed
done
