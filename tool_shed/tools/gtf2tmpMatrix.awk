#!/bin/sh
## create a initial raw count matrix from a gtf file


#Usage gtf2tmpMatrix.awk gtf($1) output($2)


awk ' BEGIN {
	FS=";"
} 
{ 
	split($1,ens,"\"")
	split($4,na,"\"")
	print ens[2] "\t" na[2]
} ' $1 | sort -u > $2