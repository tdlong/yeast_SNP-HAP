dir=$1
# get the header
head -n1 $dir/$2 >haps.txt
for fname in $dir/*.txt
do
    tail -n+2 $fname >>haps.txt
done
gzip -c haps.txt >$3

