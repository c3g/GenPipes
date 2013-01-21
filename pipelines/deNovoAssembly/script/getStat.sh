export min=$(tail -1  $1 | sed "s|\t| |g" | cut -d\  -f1)

export max=$(awk ' NR==1 {print $1} ' $1 )
export gs1=$(awk ' BEGIN {a=0;b=0} {if ($1 >= 1) {a=a+$1;b=b+1}} END { print a ";" b ";" a/b } ' $1 | cut -d\;  -f1)
export cn1=$(awk ' BEGIN {a=0;b=0} {if ($1 >= 1) {a=a+$1;b=b+1}} END { print a ";" b ";" a/b } ' $1 | cut -d\;  -f2)
export mean1=$(awk ' BEGIN {a=0;b=0} {if ($1 >= 1) {a=a+$1;b=b+1}} END { print a ";" b ";" a/b } ' $1 | cut -d\;  -f3)
#export genN1=$(awk ' {if ($1 >= 1) {split($2,nam,"_"); print nam[1] "_" nam[2] } } ' $1 | sort -u | wc -l)
export median1=$(awk -v siz=$cn1 ' BEGIN {a=0;b=siz/2;c=0} {if (c <= b) {a=a+$1;c=c+1;d=$1} else {print d ; exit} } '  $1)
export n50C1=$(awk -v siz=$gs1 ' BEGIN {a=0;b=siz/2;c=0} {if (a <= b) {a=a+$1;c=c+1;d=$1} else {print d ; exit} } '  $1)
export n50C1N=$(awk -v siz=$gs1 ' BEGIN {a=0;b=siz/2;c=0} {if (a <= b) {a=a+$1;c=c+1;d=$1} else {print c ; exit} } '  $1)
#export gen50N1=$(awk -v siz=$gs1 ' BEGIN {a=0;b=siz/2;c=0} {if (a <= b) {split($2,nam,"_"); print nam[1] "_" nam[2] } } ' $1 | sort -u | wc -l)
export n90C1=$(awk -v siz=$gs1 ' BEGIN {a=0;b=siz*0.9;c=0} {if (a <= b) {a=a+$1;c=c+1;d=$1} else {print d ; exit} } '  $1)
export n90C1N=$(awk -v siz=$gs1 ' BEGIN {a=0;b=siz*0.9;c=0} {if (a <= b) {a=a+$1;c=c+1;d=$1} else {print c ; exit} } '  $1)
#export gen90N1=$(awk -v siz=$gs1 ' BEGIN {a=0;b=siz*0.9;c=0} {if (a <= b) {split($2,nam,"_"); print nam[1] "_" nam[2] } } ' $1 | sort -u | wc -l)

export gs500=$(awk ' BEGIN {a=0;b=0} {if ($1 >= 500) {a=a+$1;b=b+1}} END { print a ";" b ";" a/b } ' $1 | cut -d\;  -f1)
export cn500=$(awk ' BEGIN {a=0;b=0} {if ($1 >= 500) {a=a+$1;b=b+1}} END { print a ";" b ";" a/b } ' $1 | cut -d\;  -f2)
#export genN500=$(awk ' {if ($1 >= 500) {split($2,nam,"_"); print nam[1] "_" nam[2] } } ' $1 | sort -u | wc -l)
export mean500=$(awk ' BEGIN {a=0;b=0} {if ($1 >= 500) {a=a+$1;b=b+1}} END { print a ";" b ";" a/b } ' $1 | cut -d\;  -f3)
export median500=$(awk -v siz=$cn500 ' BEGIN {a=0;b=siz/2;c=0} {if (c <= b) {a=a+$1;c=c+1;d=$1} else {print d ; exit} } '  $1)
export n50C500=$(awk -v siz=$gs500 ' BEGIN {a=0;b=siz/2;c=0} {if (a <= b) {a=a+$1;c=c+1;d=$1} else {print d ; exit} } '  $1)
export n50C500N=$(awk -v siz=$gs500 ' BEGIN {a=0;b=siz/2;c=0} {if (a <= b) {a=a+$1;c=c+1;d=$1} else {print c ; exit} } '  $1)
#export gen50N500=$(awk -v siz=$gs500 ' BEGIN {a=0;b=siz/2;c=0} {if (a <= b) {split($2,nam,"_"); print nam[1] "_" nam[2] } } ' $1 | sort -u | wc -l)
export n90C500=$(awk -v siz=$gs500 ' BEGIN {a=0;b=siz*0.9;c=0} {if (a <= b) {a=a+$1;c=c+1;d=$1} else {print d ; exit} } '  $1)
export n90C500N=$(awk -v siz=$gs500 ' BEGIN {a=0;b=siz*0.9;c=0} {if (a <= b) {a=a+$1;c=c+1;d=$1} else {print c ; exit} } '  $1)
#export gen90N500=$(awk -v siz=$gs500 ' BEGIN {a=0;b=siz*0.9;c=0} {if (a <= b) {split($2,nam,"_"); print nam[1] "_" nam[2] } } ' $1 | sort -u | wc -l)


export gs1000=$(awk ' BEGIN {a=0;b=0} {if ($1 >= 1000) {a=a+$1;b=b+1}} END { print a ";" b ";" a/b } ' $1 | cut -d\;  -f1)
export cn1000=$(awk ' BEGIN {a=0;b=0} {if ($1 >= 1000) {a=a+$1;b=b+1}} END { print a ";" b ";" a/b } ' $1 | cut -d\;  -f2)
#export genN1000=$(awk ' {if ($1 >= 1000) {split($2,nam,"_"); print nam[1] "_" nam[2] } } ' $1 | sort -u | wc -l)
export mean1000=$(awk ' BEGIN {a=0;b=0} {if ($1 >= 1000) {a=a+$1;b=b+1}} END { print a ";" b ";" a/b } ' $1 | cut -d\;  -f3)
export median1000=$(awk -v siz=$cn1000 ' BEGIN {a=0;b=siz/2;c=0} {if (c <= b) {a=a+$1;c=c+1;d=$1} else {print d ; exit} } '  $1)
export n50C1000=$(awk -v siz=$gs1000 ' BEGIN {a=0;b=siz/2;c=0} {if (a <= b) {a=a+$1;c=c+1;d=$1} else {print d ; exit} } '  $1)
export n50C1000N=$(awk -v siz=$gs1000 ' BEGIN {a=0;b=siz/2;c=0} {if (a <= b) {a=a+$1;c=c+1;d=$1} else {print c ; exit} } '  $1)
#export gen50N1000=$(awk -v siz=$gs1000 ' BEGIN {a=0;b=siz/2;c=0} {if (a <= b) {split($2,nam,"_"); print nam[1] "_" nam[2] } } ' $1 | sort -u | wc -l)
export n90C1000=$(awk -v siz=$gs1000 ' BEGIN {a=0;b=siz*0.9;c=0} {if (a <= b) {a=a+$1;c=c+1;d=$1} else {print d ; exit} } '  $1)
export n90C1000N=$(awk -v siz=$gs1000 ' BEGIN {a=0;b=siz*0.9;c=0} {if (a <= b) {a=a+$1;c=c+1;d=$1} else {print c ; exit} } '  $1)
#export gen90N1000=$(awk -v siz=$gs1000 ' BEGIN {a=0;b=siz*0.9;c=0} {if (a <= b) {split($2,nam,"_"); print nam[1] "_" nam[2] } } ' $1 | sort -u | wc -l)

#echo "$1;$cn1;$genN1;$gs1;$min;$mean1;$median1;$max;$n50C1;$n50C1N;$gen50N1;$n90C1;$n90C1N;$gen90N1;$cn500;$genN500;$gs500;500;$mean500;$median500;$max;$n50C500;$n50C500N;$gen50N500;$n90C500;$n90C500N;$gen90N500;$cn1000;$genN1000;$gs1000;1000;$mean1000;$median1000;$max;$n50C1000;$n50C1000N;$gen50N1000;$n90C1000;$n90C1000N;$gen90N1000"

echo "$1;$cn1;$gs1;$min;$mean1;$median1;$max;$n50C1;$n50C1N;$n90C1;$n90C1N;$cn500;$gs500;500;$mean500;$median500;$max;$n50C500;$n50C500N;$n90C500;$n90C500N;$cn1000;$gs1000;1000;$mean1000;$median1000;$max;$n50C1000;$n50C1000N;$n90C1000;$n90C1000N"
