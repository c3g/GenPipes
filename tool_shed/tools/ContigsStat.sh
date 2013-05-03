#!/bin/sh

#usage getStat.sh contigLengthfile 

export min=$(tail -1  $1 | cut -d\  -f1) 
export max=$(awk ' NR==1 ' $1 | cut -d\  -f1)  
export gs1=$(awk ' BEGIN {a=0;b=0} {if ($1 >= 1) {a=a+$1;b=b+1}} END { print a ";" b ";" a/b } ' $1 | cut -d\;  -f1) 
export cn1=$(awk ' BEGIN {a=0;b=0} {if ($1 >= 1) {a=a+$1;b=b+1}} END { print a ";" b ";" a/b } ' $1 | cut -d\;  -f2) 
export mean1=$(awk ' BEGIN {a=0;b=0} {if ($1 >= 1) {a=a+$1;b=b+1}} END { print a ";" b ";" a/b } ' $1 | cut -d\;  -f3) 
export median1=$(awk ' BEGIN {a=0;b='$cn1'/2;c=0} {if (c <= b) {a=a+$1;c=c+1;d=$1} else {print d ; exit} } '  $1) 
export n50C1=$(awk ' BEGIN {a=0;b='$gs1'/2;c=0} {if (a <= b) {a=a+$1;c=c+1;d=$1} else {print d ; exit} } '  $1) 
export n50C1N=$(awk ' BEGIN {a=0;b='$gs1'/2;c=0} {if (a <= b) {a=a+$1;c=c+1;d=$1} else {print c ; exit} } '  $1)
export n90C1=$(awk ' BEGIN {a=0;b='$gs1'*0.9;c=0} {if (a <= b) {a=a+$1;c=c+1;d=$1} else {print d ; exit} } '  $1) 
export n90C1N=$(awk ' BEGIN {a=0;b='$gs1'*0.9;c=0} {if (a <= b) {a=a+$1;c=c+1;d=$1} else {print c ; exit} } '  $1)

export gs500=$(awk ' BEGIN {a=0;b=0} {if ($1 >= 500) {a=a+$1;b=b+1}} END { print a ";" b ";" a/b } ' $1 | cut -d\;  -f1) 
export cn500=$(awk ' BEGIN {a=0;b=0} {if ($1 >= 500) {a=a+$1;b=b+1}} END { print a ";" b ";" a/b } ' $1 | cut -d\;  -f2) 

export mean500=$(awk ' BEGIN {a=0;b=0} {if ($1 >= 500) {a=a+$1;b=b+1}} END { print a ";" b ";" a/b } ' $1 | cut -d\;  -f3) 
export median500=$(awk ' BEGIN {a=0;b='$cn500'/2;c=0} {if (c <= b) {a=a+$1;c=c+1;d=$1} else {print d ; exit} } '  $1) 
export n50C500=$(awk ' BEGIN {a=0;b='$gs500'/2;c=0} {if (a <= b) {a=a+$1;c=c+1;d=$1} else {print d ; exit} } '  $1) 
export n50C500N=$(awk ' BEGIN {a=0;b='$gs500'/2;c=0} {if (a <= b) {a=a+$1;c=c+1;d=$1} else {print c ; exit} } '  $1)
export n90C500=$(awk ' BEGIN {a=0;b='$gs500'*0.9;c=0} {if (a <= b) {a=a+$1;c=c+1;d=$1} else {print d ; exit} } '  $1) 
export n90C500N=$(awk ' BEGIN {a=0;b='$gs500'*0.9;c=0} {if (a <= b) {a=a+$1;c=c+1;d=$1} else {print c ; exit} } '  $1)

export gs1000=$(awk ' BEGIN {a=0;b=0} {if ($1 >= 1000) {a=a+$1;b=b+1}} END { print a ";" b ";" a/b } ' $1 | cut -d\;  -f1) 
export cn1000=$(awk ' BEGIN {a=0;b=0} {if ($1 >= 1000) {a=a+$1;b=b+1}} END { print a ";" b ";" a/b } ' $1 | cut -d\;  -f2) 

export mean1000=$(awk ' BEGIN {a=0;b=0} {if ($1 >= 1000) {a=a+$1;b=b+1}} END { print a ";" b ";" a/b } ' $1 | cut -d\;  -f3) 
export median1000=$(awk ' BEGIN {a=0;b='$cn1000'/2;c=0} {if (c <= b) {a=a+$1;c=c+1;d=$1} else {print d ; exit} } '  $1) 
export n50C1000=$(awk ' BEGIN {a=0;b='$gs1000'/2;c=0} {if (a <= b) {a=a+$1;c=c+1;d=$1} else {print d ; exit} } '  $1) 
export n50C1000N=$(awk ' BEGIN {a=0;b='$gs1000'/2;c=0} {if (a <= b) {a=a+$1;c=c+1;d=$1} else {print c ; exit} } '  $1)
export n90C1000=$(awk ' BEGIN {a=0;b='$gs1000'*0.9;c=0} {if (a <= b) {a=a+$1;c=c+1;d=$1} else {print d ; exit} } '  $1) 
export n90C1000N=$(awk ' BEGIN {a=0;b='$gs1000'*0.9;c=0} {if (a <= b) {a=a+$1;c=c+1;d=$1} else {print c ; exit} } '  $1)


echo "Info   : cont#;Gsize;min;mean;median;max;N50;N50_cont#;N90;N90_cont#"
echo "1bp+   : $cn1;$gs1;$min;$mean1;$median1;$max;$n50C1;$n50C1N;$n90C1;$n90C1N"
echo "500bp+ : $cn500;$gs500;500;$mean500;$median500;$max;$n50C500;$n50C500N;$n90C500;$n90C500N"
echo "1kb+   : $cn1000;$gs1000;1000;$mean1000;$median1000;$max;$n50C1000;$n50C1000N;$n90C1000;$n90C1000N"