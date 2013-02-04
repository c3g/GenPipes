
#for i in `awk ' { print $1 } ' assembly/wheat/wheat.gtf ` ; do awk -v na=$i ' BEGIN {syn="unknown" } { if ($1 == na) {syn=$2}} END{ print na "\t" syn } ' blast/wheat/bestHit.txt ; done | sort -k1,1 > read_count/wheat/tmpmatrix.csv

for i in `awk ' { print $1 } ' $1 ` ; do awk -v na=$i ' BEGIN {syn="unknown" } { if ($1 == na) {syn=$2}} END{ print na "\t" syn } ' $2 ; done | sort -k1,1 > $3
