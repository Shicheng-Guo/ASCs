```
for i in `ls metadata.HC*.txt`
do
vdjtools CalcPairwiseDistances -p -m $i  $i.circos 
done

for i in `ls metadata.S*.txt`
do
vdjtools CalcPairwiseDistances -p -m $i  $i.circos 
done
```
