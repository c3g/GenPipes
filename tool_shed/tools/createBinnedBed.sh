#!/bin/sh
# 
##2011-09-20 - Mathieu Bourgey - mbourgey@genomequebec.com
#Usage createBinnedBed.sh chormosomeSize($1) binSize($2) output($3)

CHR_SIZE_FILE=$1
#Chromosome  size file : col1:chrName ; col2:size
BIN_SIZE=$2
OUTPUT_FILE=$3

awk -v binS=$BIN_SIZE ' {
	chrN=$1
	maxS=$2
	start=1
	end=binS
	while (end < maxS) {
		print chrN "\n" start "\n" end
		start=start+binS
		end=end+binS
	}
	print chrN "\n" start "\n" maxS
} ' $CHR_SIZE_FILE > $OUTPUT_FILE