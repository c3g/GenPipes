#!/bin/bash

if [ $# -ne 2 ]; then
    echo $0: usage: CreateMatA.sh myfile myRes
    exit 1
fi

input=$1
res=$2

## create file format
tail -n +2  $input > ${input}.tmp
awk -v res=${res} -F'[- \t]' '{print $1"\t"$2"\t"($2+res)}' ${input}.tmp > ${input}.bins
cut -f 2- ${input}.tmp > ${input}.tmp2
paste ${input}.bins ${input}.tmp2 > ${input}.MatA
rm -f ${input}.tmp ${input}.tmp2 ${input}.bins
