cd $1
for i in *readcount.csv ; do echo $i ; na=$(echo $i | cut -d\. -f1) ; 
grep "comp" $i | sort -k1,1 > ${na}.tmpsort ; join -1 1 -2 1 tmpmatrix.csv  ${na}.tmpsort > toto.join ; mv toto.join tmpmatrix.csv ; done; 
